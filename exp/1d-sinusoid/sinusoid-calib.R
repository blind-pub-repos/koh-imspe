#####
## Functions used to test 1d X space and 1d u space for kohIMSPE with sinusoid model
####

#########  Data generation

Comp_Model <- function(x, u){return(sin(u*(x*10)))}

bias <- function(x){return(1 - (1/3)*x - (2/3)*x^2)}

#Field Data
Field <- function(x, sd=0.06){
  y_true <- Comp_Model(x,u=pi/5)+bias(x)
  return(y_true+rnorm(length(y_true),sd=sd))
}

###### Calibration

bhat.fit <- function(Xf = NULL,yF = NULL,Ym = NULL,da = NULL,ga = NULL,clean=TRUE, eps = NULL){
  resid <- yF - Ym #residuals between field data and computer model estimates
  bhat <- newGPsep(Xf,resid,d=da$start,g=ga$start,dK=TRUE)
  if(ga$mle){
    cmle <- jmleGPsep(bhat, drange = c(da$min, da$max), grange = c(ga$min, ga$max), dab = da$ab, gab = ga$ab, maxit = 500)
  }else{
    cmle <- mleGPsep(bhat,tmin=da$min,tmax=da$max,ab=da$ab, maxit = 500)
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

calib <- function(u = NULL, Xf = NULL, yF = NULL, yMhat = NULL, fit = NULL, clean=TRUE, XFu = NULL, eps = NULL){
  
  col.replace <- which(apply(XFu, 2, function(X){any(is.na(X))}))
  u.count <- 1
  for(i in col.replace){
    XFu[,i] <- u[u.count]
    u.count <- u.count + 1
  }
  
  Ym <- predGPsep(yMhat,XFu,lite=TRUE)$mean
  cmle <- fit(Xf = Xf, yF = yF, Ym = Ym, clean=clean, eps = eps)
  return(cmle)
}

GP_Fits <- function(Xm = NULL, yM = NULL, Xf = NULL, yF = NULL, eps = sqrt(.Machine$double.eps), starts = 8, uhat = NULL, previous.params = NULL, db.max=  NULL, dmax = NULL){
  
  # Fit GPs to data
  
  if(is.null(previous.params)){ 
    d <- darg(d = list(mle= TRUE, ab = c(3/2,2)), X = Xm)
    if(is.null(db.max)){
      formals(bhat.fit)$da <- darg(d=list(mle=TRUE, ab = c(3/2,5)), X=Xf)
    }else{
      formals(bhat.fit)$da <- darg(d=list(mle=TRUE, ab = c(3/2,5), max = db.max), X=Xf)
    }
    formals(bhat.fit)$ga <- garg(g=list(mle=TRUE, ab = c(3/2, 7)), y=yF)
  }else{
    d <- darg(d = list(mle= TRUE, ab = c(3/2,2)), X = Xm)
    d$start <- previous.params$theta_m
    if(is.null(db.max)){
      temp <- darg(d=list(mle=TRUE, ab = c(3/2,5)), X=Xf)
    } else{
      temp <- darg(d=list(mle=TRUE, ab = c(3/2,5), max = db.max), X=Xf)
    }
    temp$start <- previous.params$theta_b
    formals(bhat.fit)$da  <- temp
    formals(bhat.fit)$ga <- garg(g=list(mle=TRUE, ab = c(3/2, 7), start = previous.params$gb), y=yF)
  }
  
  if(!is.null(dmax)){
    d$max <- dmax
  }
  
  yMhat <- newGPsep(Xm, yM, d=d$start, g = eps, dK=TRUE) #Fix nugget b/c deterministic model
  mle <- mleGPsep(yMhat, param = 'd', tmin=eps, tmax =  d$max, ab = d$ab, maxit = 500) #Estimate the lengthscale parameters
  
  
  lprior <- function(u, shape1=2, shape2=2, lwr=0, upr=1){
    sum(dbeta(u, shape1, shape2, log=TRUE))
  }
  
  
  #Optimize u
  obj <- function(x = NULL, Xf = NULL, yF = NULL, yMhat = NULL, fit = NULL, XFu = NULL, eps = NULL) {
    if(any(x < 0) || any(x > 1)) return(Inf)
    calib(u = x, Xf = Xf,yF = yF, yMhat = yMhat, fit = fit, XFu = XFu, eps = eps)$nll - lprior(x)
  }
  
  uhats <- data.frame(matrix(NA, ncol = ncol(Xm)-ncol(Xf) + 1, nrow = starts))
  names(uhats) <- c("value", paste("uhat", 1:(ncol(Xm)-ncol(Xf)), sep = ""))
  
  proposed.uhats <- tgp::lhs(starts, rect = matrix(c(0,1), ncol = 2, nrow = ncol(Xm)-ncol(Xf), byrow = TRUE))
  
  XFu <- cbind(Xf, matrix(NA,nrow = nrow(Xf), ncol=ncol(Xm)-ncol(Xf)))
  
  
  if(is.null(uhat)){
    
    for(i in 1:starts){
      soln <-  optim(proposed.uhats[i,],obj, method = "Brent", lower = 0, upper = 1,Xf=Xf, yF = yF, yMhat = yMhat, fit=bhat.fit, XFu = XFu, eps = eps)
      uhats[i,] <- c(soln$value, soln$par)
    }
    
    uhat <-  uhats[which.min(uhats$value),]
    uhat$value <- NULL
    uhat <- unlist(unname(uhat))
  }
  
  XFu <- cbind(Xf, matrix(rep(uhat,nrow(Xf)), ncol=length(uhat),byrow = TRUE))
  
  #Obtain GP estimate and predictions for bias model
  bhat <- calib(u = uhat,Xf = Xf,yF = yF,yMhat = yMhat,fit = bhat.fit, XFu = XFu ,clean=FALSE) #Returns the GP
  
  
  # Extract Scale Estimates
  KM <- covar.sep(Xm, d=mle$d, g=eps) #Correlation matrix of Comp Design Points
  tau2M <- drop(crossprod(yM, solve(KM, yM)))/length(yM)
  
  #Covariance for Bias Model at Field Locations
  #Correlation Matrix
  
  LS.B <- unlist(unname(bhat[,1:ncol(Xf)]))
  Nug.B <- unlist(unname(bhat[,ncol(Xf) + 1]))
  
  KB <- covar.sep(Xf, d=LS.B, g=Nug.B)
  #Pred of comp model at field locations
  XFuhat <- cbind(Xf,matrix(rep(uhat,nrow(Xf)),ncol=length(uhat),byrow=TRUE))
  Ym <- predGPsep(yMhat,XFuhat,lite=TRUE)$mean
  #Difference between Field obs and Comp model estimates at field locations
  Yb <- yF - Ym
  #Estimate Tau2 for Bias Model
  tau2B <- drop(crossprod(Yb,solve(KB,Yb))) / length(Yb)
  
  #Store Results
  LS.M <- mle$d
  
  deleteGPsep(yMhat)
  deleteGPsep(bhat$gp)
  
  return(list(uhat=uhat, tau2M=tau2M, LS.M = LS.M, tau2B=tau2B, LS.B = LS.B, Nug.B = Nug.B))
}


