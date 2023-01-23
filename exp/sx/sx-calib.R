### Calibration Example: SX Application

Comp_Model <- function(X,U, ranges, t = NULL, h= NULL, outs = 1:3, log = FALSE){
  
  y.out <- RE.conc(X = X, U = U, ranges = ranges, t = t, h = h)
  if(is.vector(y.out)){
    if(log == TRUE){
      y.out <- log(y.out[outs])
    }else{
      y.out <- y.out[outs]
    }
  }else{
    if(log == TRUE){
      y.out <- drop(log(y.out[, outs]))
    }else{
      y.out <- drop(y.out[, outs])
    }
  }
  return(y.out)
}


bhat.fit.sep <- function(Xf = NULL,yF = NULL,Ym = NULL,da = NULL,ga = NULL,clean=TRUE, eps = NULL){
  resid <- yF - Ym #residuals between field data and computer model estimates
  
  bhat <- newGPsep(Xf,resid,d=da$start,g=ga$start,dK=TRUE)
  if(ga$mle){
    cmle <- jmleGPsep(bhat, drange = c(da$min, da$max), grange = c(ga$min, ga$max), dab = da$ab, gab = ga$ab, maxit = 500)
  }else{
    cmle <- jmleGPsep(bhat,tmin=da$min,tmax=da$max,ab=da$ab, maxit = 500)
  }
  cmle$nll <- - llikGPsep(bhat, da=da$ab,gab=ga$ab)
  if(clean){
    deleteGPsep(bhat)
  }else{
    cmle$gp <- bhat
    cmle$gptype <- 'sep'
  }
  return(cmle)
}

bhat.fit.iso <- function(Xf = NULL,yF = NULL,Ym = NULL,da = NULL,ga = NULL,clean=TRUE, eps = NULL){
  resid <- yF - Ym #residuals between field data and computer model estimates
  
  bhat <- newGP(Xf,resid,d=da$start,g=ga$start,dK=TRUE)
  if(ga$mle){
    cmle <- jmleGP(bhat, drange = c(da$min, da$max), grange = c(ga$min, ga$max), dab = da$ab, gab = ga$ab)
  }else{
    cmle <- jmleGP(bhat,tmin=da$min,tmax=da$max,ab=da$ab)
  }
  cmle$nll <- - llikGP(bhat, da=da$ab,gab=ga$ab)
  if(clean){
    deleteGP(bhat)
  }else{
    cmle$gp <- bhat
    cmle$gptype <- 'iso'
  }
  return(cmle)
}

calib <- function(u = NULL, Xf = NULL, yF = NULL, Xm  = NULL, yMhat = NULL, yMhat.mle = NULL, fit = NULL, clean=TRUE, XFu = NULL, eps = NULL, KiyM = NULL){
  
  for(i in (ncol(Xf)+1):ncol(XFu)){
    XFu[,i] <- u[i-ncol(Xf)]
  }
  
  Ym <- drop(cov_gen(XFu, Xm, theta = yMhat.mle$d, type = "Gaussian") %*% KiyM)
  
  #Ym <- predGPsep(yMhat,XFu,lite=TRUE)$mean
  cmle <- fit(Xf = Xf, yF = yF, Ym = Ym, clean=clean, eps = eps)
  return(cmle)
}

GP_Fits <- function(Xm = NULL, yM = NULL, Xf = NULL, yF = NULL, eps = sqrt(.Machine$double.eps), starts = 10, uhat = NULL, previous.params = NULL, sep = NULL){
  
  # Fit GPs to data
  
  if(sep == TRUE){
    bhat.fit <- bhat.fit.sep
    
    if(is.null(previous.params)){ 
      d <- darg(d = list(mle= TRUE, ab = c(3/2,0.9), max = 10), X = Xm)
      formals(bhat.fit)$da <- darg(d=list(mle=TRUE, ab = c(3/2,5/2), max = 0.25), X=Xf)
      formals(bhat.fit)$ga <- garg(g=list(mle=TRUE, ab = c(3/2, 0.05), max = 2*var(yF)), y=yF)
    }else{
      d <- darg(d = list(mle= TRUE, ab = c(3/2,0.9), max = 10), X = Xm)
      d$start <- previous.params$theta_m
      temp <- darg(d=list(mle=TRUE, ab = c(3/2,0.9), max = 0.25), X=Xf)
      temp$start <- previous.params$theta_b
      formals(bhat.fit)$da  <- temp
      formals(bhat.fit)$ga <- garg(g=list(mle=TRUE, ab = c(3/2, 0.05), start = previous.params$gb, max = 2*var(yF)), y=yF)
    }
    
  }else if(sep == FALSE){
    bhat.fit <- bhat.fit.iso
    
    if(is.null(previous.params)){ 
      d <- darg(d = list(mle= TRUE, ab = c(3/2,5/4), max = 10), X = Xm)
      formals(bhat.fit)$da <- darg(d=list(mle=TRUE, ab = c(3/2,5/2), max = 0.25), X=Xf)
      formals(bhat.fit)$ga <- garg(g=list(mle=TRUE, ab = c(3/2, 0.05), max = 2*var(yF)), y=yF)
    }else{
      d <- darg(d = list(mle= TRUE, ab = c(3/2,5/4), max = 10), X = Xm)
      d$start <- previous.params$theta_m
      temp <- darg(d=list(mle=TRUE, ab = c(3/2,5/2), max = 0.25), X=Xf)
      temp$start <- previous.params$theta_b
      formals(bhat.fit)$da  <- temp
      formals(bhat.fit)$ga <- garg(g=list(mle=TRUE, ab = c(3/2, 0.05), start = previous.params$gb, max = 2*var(yF)), y=yF)
    }
    
  }
  
  
  yMhat <- newGPsep(Xm, yM, d=d$start, g = eps, dK=TRUE) #Fix nugget b/c deterministic model
  mle <- mleGPsep(yMhat, param = 'd', tmin=eps, tmax =  5, ab = d$ab, maxit = 1000) #Estimate the lengthscale parameters
  
  Km <- cov_gen(Xm, theta = mle$d, type = "Gaussian")
  diag(Km) <- diag(Km) + eps
  KiyM <- solve(Km, yM)
  
  lprior <- function(u, shape1=2, shape2=2) 
  {
    sum(dbeta(u, shape1, shape2, log=TRUE))
  }
  
  
  #Optimize u
  #obj <- function(x = NULL, Xf = NULL, yF = NULL, yMhat = NULL, fit = NULL, XFu = NULL, eps = NULL) {
  obj <- function(x = NULL, Xf = NULL, yF = NULL, Xm = NULL, fit = NULL, yMhat.mle = NULL, XFu = NULL, eps = NULL, KiyM = NULL) {
    if(any(x < 0) || any(x > 1)) return(Inf)
    calib(u = x, Xf = Xf,yF = yF, Xm = Xm, yMhat.mle = yMhat.mle, fit = fit, XFu = XFu, eps = eps, KiyM = KiyM)$nll - lprior(x)
  }
  
  
  uhats <- data.frame(matrix(NA, ncol = ncol(Xm)-ncol(Xf) + 1, nrow = starts))
  names(uhats) <- c("value", paste("uhat", 1:(ncol(Xm)-ncol(Xf)), sep = ""))
  
  proposed.uhats <- randomLHS(n = starts, k = ncol(Xm)-ncol(Xf))
  XFu <- cbind(Xf, matrix(NA,nrow = nrow(Xf), ncol=ncol(Xm)-ncol(Xf)))
  
  if(is.null(uhat)){
    
    for(i in 1:starts){
      soln <-  optim(proposed.uhats[i,],obj, method = "Nelder-Mead", Xf=Xf, yF = yF, Xm = Xm, fit = bhat.fit, yMhat.mle = mle, XFu = XFu, eps = eps, KiyM = KiyM)
      uhats[i,] <- c(soln$value, soln$par)
    }
    uhat <-  uhats[which.min(uhats$value),]
    uhat$value <- NULL
    uhat <- unlist(unname(uhat))
  }
  
  XFu <- cbind(Xf, matrix(rep(uhat,nrow(Xf)), ncol=length(uhat),byrow = TRUE))
  
  #Obtain GP estimate and predictions for bias model
  bhat <- calib(u = uhat,Xf = Xf,yF = yF, Xm = Xm, yMhat.mle = mle,fit = bhat.fit, XFu = XFu ,clean=FALSE, KiyM = KiyM) #Returns the GP
  
  
  # Extract Scale Estimates
  KM <- covar.sep(Xm, d=mle$d, g=eps) #Correlation matrix of Comp Design Points
  tau2M <- drop(crossprod(yM, solve(KM, yM)))/length(yM)
  
  #Covariance for Bias Model at Field Locations
  #Correlation Matrix
  
  if(sep == TRUE){
    LS.B <- unlist(unname(bhat[,1:ncol(Xf)]))
    Nug.B <- unlist(unname(bhat[,ncol(Xf) + 1]))
    deleteGPsep(bhat$gp)
  }else if(sep == FALSE){
    LS.B <- rep(bhat$d, times = ncol(Xf))
    Nug.B <- bhat$g
    deleteGP(bhat$gp)
  }
  
  KB <- covar.sep(Xf, d=LS.B, g=Nug.B)
  #Pred of comp model at field locations
  XFuhat <- cbind(Xf,matrix(rep(uhat,nrow(Xf)),ncol=length(uhat),byrow=TRUE))
  
  Ym <- drop(cov_gen(XFuhat, Xm, theta = mle$d, type = "Gaussian") %*% KiyM)
  
  #Ym <- predGPsep(yMhat,XFuhat,lite=TRUE)$mean
  #Difference between Field obs and Comp model estimates at field locations
  Yb <- yF - Ym
  #Estimate Tau2 for Bias Model
  tau2B <- drop(crossprod(Yb,solve(KB,Yb))) / length(Yb)
  
  #Store Results
  LS.M <- mle$d
  
  deleteGPsep(yMhat)

  return(list(uhat=uhat, tau2M=tau2M, LS.M = LS.M, tau2B=tau2B, LS.B = LS.B, Nug.B = Nug.B))
}
