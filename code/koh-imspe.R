
############ IMSPE Functions  #################

### Use this function to find minimizing IMSPE model point.  
### After modular KOH calibration assuming a deterministic simulator, Gaussian covariance functions, and homoskedastic field noise, supply the following vars:

## Xm  Design for observed simulator data (matrix)
## Xf  Design for observed field data (matrix)
## uhat MLE estimation of simulator calibration parameter (may only currently support a 1-D uhat...not sure, but it would be easy to fix)
## Vb  Scale hyperparameter for the bias GP
## Vm  Scale hyperparameter for the model GP
## theta_m  lengthscale hyperparameters for the model GP
## theta_b  lengthscale hyperparameters for the bias GP
## gb  nugget for the gias GP
## type  Covariance function for GP fits (note, currently only "Gaussian" is supported)
## p  Used for cokriging, number of columns of the response (currently only support for p = 1)
## starts  number of random starts of optim in the search for minimum IMSPE
## eps  value added to diagonal of matricies before inversion.

IMSPE.optim <- function(Xm = NULL, Xf = NULL, uhat = NULL, Vb = NULL, Vm = NULL, theta_m = NULL, theta_b = NULL,
                        gb = NULL, type = NULL, starts = 100,candidate = TRUE, eps = sqrt(.Machine$double.eps), log.IMSPE = TRUE, optim.control = list(pgtol = 0.05, trace = 0, lmm = 13)){
  
  nm <- nrow(Xm)
  nf <- nrow(Xf)

  
  Xfm <- rbind(cbind(Xf,rep(1, nrow(Xf)) %x% t(uhat)), Xm)
  Xfuhat <- Xfm[1:nf,]
  Km <- cov_gen(Xfm, theta = theta_m, type = type)
  Km <- Vm*(Km + diag(nrow(Km))*eps)
  
  Km[1:nf,1:nf] <- Km[1:nf,1:nf] + (Vb)*(cov_gen(Xf, theta = theta_b, type = type) + diag(nf)* gb)
  
  Ki <- solve(Km)
  
  Ki <- (Ki + t(Ki))/2
  
  if(type == "Gaussian"){
    W11.fn <- W11.gaussian
    W12.fn <- W12.gaussian
    W11.xtilde <- W11.xtilde.gaussian
    W12.xtilde <- W12.xtilde.gaussian
  }
  
  W11_mr1 <- W11.fn(Xf = Xf, uhat = uhat, Xm = Xm, theta_m=  theta_m, type = type)
  W12 <- W12.fn(Xf = Xf, uhat = uhat, Xm = Xm, theta_m=  theta_m, theta_b = theta_b, Xfuhat = Xfuhat, type = type)
  W22 <- matrix(0, nrow = nrow(Xf) + nrow(Xm) + 1, ncol = nrow(Xf) + nrow(Xm) + 1)
  W22[1:nrow(Xf), 1:nrow(Xf)] <- hetGP::Wij(mu1 = Xf, theta = theta_b, type =type)
  
  params <- list()
  params$Xm <- Xm
  params$Xf <- Xf
  params$Ki <- Ki
  params$W11_mr1 <- W11_mr1
  params$W12 <- W12
  params$W22 <- W22
  params$type <- type
  params$eps <- eps
  params$Vm <- Vm
  params$Vb <- Vb
  params$theta_m <- theta_m
  params$theta_b <- theta_b
  params$uhat <- uhat
  params$Km <- Km
  params$log.IMSPE <- log.IMSPE
  params$Xfuhat <- Xfuhat
  params$Xfm <- Xfm
  
  params.cand <- params
  params.cand$log.IMSPE <- FALSE
  
  inits <- randomLHS(n = starts, k = ncol(Xm))
  imspes <- apply(inits,1,IMSPE.obj, params = params.cand)
  
  if(is.null(optim.control$pgtol)){
    optim.control$pgtol <- 0.1
  }
  
  # if(optim.control$pgtol == "adaptive"){
  #   dimspes <- abs(apply(inits,1,dIMSPE, params = params))
  #   optim.control$pgtol <- unname(max(abs(quantile(as.vector(dimspes), probs = c(0.025, 0.975)))))
  #   optim.control.save <<- c(optim.control$pgtol, optim.control.save)
  # }
  
  if(is.null(optim.control$trace)){
    optim.control$trace <- 0
  }
  
  if(is.null(optim.control$lmm)){
    optim.control$lmm <- 13
  }

  
  if(any(imspes < 0)){
    warning("negative IMSPE calculated in candidate set", .immediate = TRUE)
  }
  
  if(candidate){
    
    x.tilde <- inits[which.min(imspes),]
    
  }else{
  
    cands <- sort(imspes, index.return = TRUE)$ix
    inits <- inits[cands[1:ceiling(length(imspes)*0.075)],]
    imspes <- imspes[cands[1:ceiling(length(imspes)*0.075)]]
    
  candidates <- data.frame(matrix(NA, ncol = ncol(Xm)+1, nrow = nrow(inits)))
  names(candidates) <- c("value", paste("par", 1:ncol(Xm), sep = ""))
  
  for(i in 1:nrow(inits)){
    
    out <- optim(par = inits[i,], fn = IMSPE.obj, gr = dIMSPE, method = "L-BFGS-B", 
                 lower = rep(0, times = ncol(Xm)), upper = rep(1, times = ncol(Xm)), 
                 control = optim.control, params = params)

    if(out$convergence != 0){
      warning(paste("Optim for IMSPE failed to converge:  Code", out$convergence, "increase pgtol value", sep = " ", collapse = " "), immediate. = TRUE)
      if(imspes[i] <= 0){
        candidates[i,] <- NA
      }else{
        if(log.IMSPE){
        candidates[i,] <- c(log(imspes[i]), inits[i,])
        }else{
          candidates[i,] <- c(imspes[i], inits[i,])
        }
      }
    }else{
      if(log.IMSPE == FALSE & out$value < 0){
        candidates[i, ] <- NA
      }else{
        candidates[i,] <- c(out$value, out$par)
      }
     }
  }
  
  x.tilde <- unlist(unname(drop(candidates[which.min(candidates$value),-1])))
  }
  
  if(length(x.tilde) == 0){
    stop("No x.tilde found, try more random starts but likely numerical issues with IMSPE calculation.")
  }
  
  return(x.tilde)
}

######################################################
## Gaussian Kernel Functions ############################
##############################################################

W12.gaussian <- function(Xf = NULL, uhat = NULL, Xm = NULL, theta_m = NULL, theta_b = NULL, Xfuhat = NULL, type = "Gaussian"){
  
  ## Iterate over model fit data
  
  W.temp <- matrix(0, nrow = (nrow(Xf) + nrow(Xm)), ncol = nrow(Xf) + nrow(Xm))
  
  theta.i <- theta_b
  theta.j <- theta_m
  
  for(i in 1:(nrow(Xf))){
    # Iterate over bias fit data
    W.temp[i,] <- 1
    xi <- Xf[i,]
    for(j in 1:(nrow(Xf)+ nrow(Xm))){
      ## xi is from the bias kernel and xj is from the model kernel
      
      if(j <= nrow(Xf)){
        xj <- Xfuhat[j,]
      }else{
        xj <- Xm[j-nrow(Xf),]
      }
      for(l in 1:length(theta_b)){
        arg1 <- (theta.j[l]*xi[l] + theta.i[l]*xj[l])/(theta.i[l] + theta.j[l])
        arg2 <- (arg1-1)/sqrt(theta.i[l]*theta.j[l]/(theta.i[l] + theta.j[l]))
        arg1 <- arg1/sqrt(theta.i[l]*theta.j[l]/(theta.i[l] + theta.j[l]))
        arg <- erf(arg1) - erf(arg2)
        
        arg <- exp(-(xi[l]-xj[l])^2/(theta.i[l] + theta.j[l]))*(0.5*sqrt(pi * (theta.i[l]*theta.j[l]/(theta.i[l]+theta.j[l]))))*arg
        
        W.temp[i,j] <- W.temp[i,j]*arg
      }
      
      if(j > nrow(Xf)){
        for(l in (length(theta_b)+1):length(theta_m)){
          arg <- exp(-(uhat[l-length(theta_b)]- xj[l])^2/theta.j[l])
          W.temp[i,j] <- W.temp[i,j]*arg
        }
      }
    }
  }
  
  
  return(W.temp)
  
}

W12.xtilde.gaussian <- function(W12 = NULL, xtilde = NULL, theta_b = NULL, theta_m = NULL, uhat = NULL, Xf = NULL, Xm = NULL, type = "Gaussian"){
  
  W.temp <- rep(1, times = nrow(Xf))
  
  theta.i <- theta_b
  theta.j <- theta_m
  u.addtl <- 0
  
  for(l in (length(theta_b)+1):length(theta_m)){
    u.addtl <- -(uhat[l-length(theta_b)]- xtilde[l])^2/theta.j[l] + u.addtl
  }
  u.addtl <- exp(u.addtl)
  
  for(i in 1:(nrow(Xf))){
    # Iterate over bias fit data
    xi <- Xf[i,]
    for(l in 1:length(theta_b)){
      arg1 <- (theta.j[l]*xi[l] + theta.i[l]*xtilde[l])/(theta.i[l] + theta.j[l])
      arg2 <- (arg1-1)/sqrt(theta.i[l]*theta.j[l]/(theta.i[l] + theta.j[l]))
      arg1 <- arg1/sqrt(theta.i[l]*theta.j[l]/(theta.i[l] + theta.j[l]))
      arg <- erf(arg1) - erf(arg2)
      
      arg <- exp(-(xi[l]-xtilde[l])^2/(theta.i[l] + theta.j[l]))*(0.5*sqrt(pi * (theta.i[l]*theta.j[l]/(theta.i[l]+theta.j[l]))))*arg
      
      W.temp[i] <- W.temp[i]*arg
    }
  }
  
  W.temp <- W.temp*u.addtl
  
  return(rbind(cbind(W12, c((W.temp), rep(0, times = nrow(Xm)))), rep(0, times = nrow(Xf) + nrow(Xm)+1)))
  
}

W11.gaussian <- function(Xf = NULL, uhat = NULL, Xm = NULL, theta_m =  NULL, type = "Gaussian"){
  
  ## Iterate over model fit data
  
  Xfmb <- rbind(Xf, Xm[,1:ncol(Xf), drop = FALSE])
  
  ## all elements related to X space
  
  W.temp <- hetGP::Wij(Xfmb, theta = theta_m[1:ncol(Xf)], type = type)
  
  ## exponentials 
  Wu.mm <- cov_gen(Xm[,(ncol(Xf)+1):ncol(Xm), drop = FALSE], matrix(uhat, ncol = length(uhat), nrow = 1), 
                     theta = theta_m[(ncol(Xf)+1):ncol(Xm)])

  Wu.MM <- tcrossprod(Wu.mm)
  
  ## elements of model data only, adding values related to u space
  W.temp[(nrow(Xf) + 1):nrow(W.temp),(nrow(Xf) + 1):nrow(W.temp)] <- W.temp[(nrow(Xf) + 1):nrow(W.temp),(nrow(Xf) + 1):nrow(W.temp)] *Wu.MM
  
  ## Elements related to kernel distances between field data and model data, adding values related to u space
  W.temp[1:nrow(Xf),(nrow(Xf)+1):nrow(W.temp)] <-  W.temp[1:nrow(Xf),(nrow(Xf)+1):nrow(W.temp)] * matrix(Wu.mm, nrow = nrow(Xf), ncol = nrow(Xm), byrow = TRUE)
  W.temp[(nrow(Xf) + 1):nrow(W.temp), 1:nrow(Xf)] <- t(W.temp[1:nrow(Xf),(nrow(Xf)+1):nrow(W.temp)])
  
  ## Note that because field data and predictive location are both located at uhat, covariance is 1 at elements only related to field data
  
  return(W.temp)
}

W11.xtilde.gaussian <- function(xtilde = NULL, W11_mr1 = NULL, Xf = NULL, Xm = NULL, uhat=  NULL, theta_m = NULL, type = "Gaussian"){
  
  Xfmb <- rbind(Xf, Xm[,1:ncol(Xf), drop = FALSE])
  W.temp <- hetGP::Wij(matrix(xtilde[1:ncol(Xf)], nrow = 1, ncol= ncol(Xf)), mu2 = Xfmb, theta = theta_m[1:ncol(Xf)], type = type)
  uhatutilde <- exp(- sum(((uhat-xtilde[(ncol(Xf)+1):length(xtilde)])^2)/theta_m[(ncol(Xf)+1):length(xtilde)]))
  
  W.temp <- W.temp * uhatutilde
  
  ### Integral related to multiplying kernels related to xtilde and model data
  W.temp[1,(nrow(Xf)+1):ncol(W.temp)] <- W.temp[1,(nrow(Xf)+1):ncol(W.temp), drop=  FALSE] * hetGP::cov_gen(matrix(uhat, nrow = 1, ncol= length(uhat)),Xm[,(ncol(Xf)+1):ncol(Xm), drop= FALSE], theta =theta_m[(ncol(Xf)+1):length(theta_m)], type = type)
  
  
  w.temp <- hetGP::Wij(matrix(xtilde[1:ncol(Xf)], ncol = ncol(Xf), nrow=  1), theta = theta_m[1:ncol(Xf)], type = type) * uhatutilde^2
  
  return(cbind(rbind(W11_mr1, W.temp), rbind(t(W.temp), w.temp)))
  
}


W12.xtilde <- function(W12 = NULL, xtilde = NULL, theta_b = NULL, theta_m = NULL, uhat = NULL, Xf = NULL, Xm = NULL, type = NULL){
  
  W.temp <- rep(1, times = nrow(Xf))
  
  theta.i <- theta_b
  theta.j <- theta_m
  u.addtl <- 0
  
  for(l in (length(theta_b)+1):length(theta_m)){
    u.addtl <- -(uhat[l-length(theta_b)]- xtilde[l])^2/theta.j[l] + u.addtl
  }
  u.addtl <- exp(u.addtl)
  
  for(i in 1:(nrow(Xf))){
    # Iterate over bias fit data
    xi <- Xf[i,]
    for(l in 1:length(theta_b)){
      arg1 <- (theta.j[l]*xi[l] + theta.i[l]*xtilde[l])/(theta.i[l] + theta.j[l])
      arg2 <- (arg1-1)/sqrt(theta.i[l]*theta.j[l]/(theta.i[l] + theta.j[l]))
      arg1 <- arg1/sqrt(theta.i[l]*theta.j[l]/(theta.i[l] + theta.j[l]))
      arg <- erf(arg1) - erf(arg2)
      
      arg <- exp(-(xi[l]-xtilde[l])^2/(theta.i[l] + theta.j[l]))*(0.5*sqrt(pi * (theta.i[l]*theta.j[l]/(theta.i[l]+theta.j[l]))))*arg
      
      W.temp[i] <- W.temp[i]*arg
    }
  }
  
  W.temp <- W.temp*u.addtl
  
  return(rbind(cbind(W12, c((W.temp), rep(0, times = nrow(Xm)))), rep(0, times = nrow(Xf) + nrow(Xm)+1)))
  
}

W11_mr1.fn <- function(Xf = NULL, uhat = NULL, Xm = NULL, theta_m =  NULL, type = NULL){
  
  ## Iterate over model fit data
  
  Xfmb <- rbind(Xf, Xm[,1:ncol(Xf), drop = FALSE])
  
  ## all elements related to X space
  
  W.temp <- hetGP::Wij(Xfmb, theta = theta_m[1:ncol(Xf)], type = type)
  
  ## exponentials 
  Wu.mm <- cov_gen(Xm[,(ncol(Xf)+1):ncol(Xm), drop = FALSE], matrix(uhat, ncol = length(uhat), nrow = 1), 
                   theta = theta_m[(ncol(Xf)+1):ncol(Xm)])
  
  Wu.MM <- tcrossprod(Wu.mm)
  
  ## elements of model data only, adding values related to u space
  W.temp[(nrow(Xf) + 1):nrow(W.temp),(nrow(Xf) + 1):nrow(W.temp)] <- W.temp[(nrow(Xf) + 1):nrow(W.temp),(nrow(Xf) + 1):nrow(W.temp)] *Wu.MM
  
  ## Elements related to kernel distances between field data and model data, adding values related to u space
  W.temp[1:nrow(Xf),(nrow(Xf)+1):nrow(W.temp)] <-  W.temp[1:nrow(Xf),(nrow(Xf)+1):nrow(W.temp)] * matrix(Wu.mm, nrow = nrow(Xf), ncol = nrow(Xm), byrow = TRUE)
  W.temp[(nrow(Xf) + 1):nrow(W.temp), 1:nrow(Xf)] <- t(W.temp[1:nrow(Xf),(nrow(Xf)+1):nrow(W.temp)])
  
  ## Note that because field data and predictive location are both located at uhat, covariance is 1 at elements only related to field data
  
  return(W.temp)
}

W11.xtilde <- function(xtilde = NULL, W11_mr1 = NULL, Xf = NULL, Xm = NULL, uhat=  NULL, theta_m = NULL, type = NULL){
  
  Xfmb <- rbind(Xf, Xm[,1:ncol(Xf), drop = FALSE])
  W.temp <- hetGP::Wij(matrix(xtilde[1:ncol(Xf)], nrow = 1, ncol= ncol(Xf)), mu2 = Xfmb, theta = theta_m[1:ncol(Xf)], type = type)
  uhatutilde <- exp(- sum(((uhat-xtilde[(ncol(Xf)+1):length(xtilde)])^2)/theta_m[(ncol(Xf)+1):length(xtilde)]))
  
  W.temp <- W.temp * uhatutilde
  
  ### Integral related to multiplying kernels related to xtilde and model data
  W.temp[1,(nrow(Xf)+1):ncol(W.temp)] <- W.temp[1,(nrow(Xf)+1):ncol(W.temp), drop=  FALSE] * hetGP::cov_gen(matrix(uhat, nrow = 1, ncol= length(uhat)),Xm[,(ncol(Xf)+1):ncol(Xm), drop= FALSE], theta =theta_m[(ncol(Xf)+1):length(theta_m)], type = type)
  
  
  w.temp <- hetGP::Wij(matrix(xtilde[1:ncol(Xf)], ncol = ncol(Xf), nrow=  1), theta = theta_m[1:ncol(Xf)], type = type) * uhatutilde^2
  
  return(cbind(rbind(W11_mr1, W.temp), rbind(t(W.temp), w.temp)))
  
}





#########################################################################
###########################################################################
################################################################################3

Kn1i.fn <- function(Ki = NULL, b = NULL, K_x.fm = NULL, Kikx = NULL){
  
  ## top left
  
  Kn1i <- matrix(NA, ncol = ncol(Ki) +1, nrow = nrow(Ki)+1)
  
  tl <- Ki + tcrossprod(Kikx) /b
  Kn1i[1:nrow(Ki), 1:nrow(Ki)] <- tl
  

  
  ## top right
  
  tr <- -Kikx/b
  Kn1i[1:nrow(Ki), nrow(Ki)+1] <- tr
  Kn1i[nrow(Ki)+1, 1:nrow(Ki)] <- t(tr)
  Kn1i[nrow(Ki)+1, nrow(Ki)+1] <- 1/b
  
  return(Kn1i)
  
}

IMSPE.obj <- function(par,params = NULL){
  
  xtilde <- matrix(par, nrow = 1, ncol= length(par))
  
  
  Xf <- params$Xf
  Xm <- params$Xm
  uhat <- params$uhat
  theta_m <- params$theta_m
  Vm <- params$Vm
  Vb <- params$Vb
  Ki <- params$Ki
  W12 <- params$W12
  W22 <- params$W22
  theta_b <- params$theta_b
  W11_mr1 <- params$W11_mr1
  type <- params$type
  eps <- params$eps
  Km <- params$Km
  log.IMSPE <- params$log.IMSPE
  Xfuhat <- params$Xfuhat
  Xfm <- params$Xfm
  
  
  K_x.fm <- Vm * cov_gen(xtilde, Xfm, theta= theta_m, type = type)
  
  Kikx <- solve(Km,t(K_x.fm))
  
  
  b <-  drop(Vm*(1 + eps) - K_x.fm %*% Kikx)
  
  Kn1i <- Kn1i.fn(Ki = Ki, b = b, K_x.fm = K_x.fm, Kikx = Kikx)
  
  
  
  ## Update all W matricies
  
  W12.x <- W12.xtilde(W12 = W12, xtilde = par, theta_b = theta_b, theta_m = theta_m, uhat = uhat, Xf = Xf, Xm = Xm, type = type)
  W11.x <- W11.xtilde(xtilde = par, W11_mr1 = W11_mr1, Xf = Xf, Xm = Xm, uhat=  uhat, theta_m = theta_m, type = type)
  
  temp <- sum(Kn1i * (W11.x * Vm^2)) + 2*sum(Kn1i *(W12.x*Vm * Vb)) + sum(Kn1i*(W22 * Vb^2))
  
    temp <- Vb + Vm*(1+eps) - temp
    
    if(log.IMSPE){
      if(temp <= 0){
        stop("Negative values for IMSPE encountered, unable to calculate log(IMSPE).  Try using large candidate based search.  Possibly little improvement can be gained from additional model points")
      }else{
        temp <- log(temp)
        }
    }
  
  return(temp)
  
}


############ Gradient Functions ################

dIMSPE <- function(par, params = NULL){
  
  xtilde <- par
  
  Xf <- params$Xf
  Xm <- params$Xm
  uhat <- params$uhat
  theta_m <- params$theta_m
  Vm <- params$Vm
  Vb <- params$Vb
  Ki <- params$Ki
  Km <- params$Km
  W12 <- params$W12
  W22 <- params$W22
  theta_b <- params$theta_b
  W11_mr1 <- params$W11_mr1
  type <- params$type
  eps <- params$eps
  log.IMSPE <- params$log.IMSPE
  Xfuhat <- params$Xfuhat
  Xfm <- params$Xfm
 
  K_x.fm <- cov_gen(matrix(xtilde, nrow=  1, ncol= length(xtilde)), Xfm, theta= theta_m, type = type)*Vm
  
  
  Kikx <- solve(Km,t(K_x.fm))
  
  
  b <-  drop(Vm*(1 + eps) - K_x.fm %*% Kikx)
  
  Kn1i <- Kn1i.fn(Ki = Ki, b = b, K_x.fm = K_x.fm, Kikx = Kikx)
  
  
  dKi.dxtilde <- dKi(xtilde = xtilde, uhat = uhat, Xm = Xm, Xf = Xf, Xfuhat = Xfuhat, theta_m = theta_m, Ki = Ki,  Vm = Vm, Km = Km, Kikx = Kikx ,K_x.fm= K_x.fm, b = b,
                     eps = eps, type = type)
  
  ## Update all W matricies
  
  W12.x <- W12.xtilde(W12 = W12, xtilde = par, theta_b = theta_b, theta_m = theta_m, uhat = uhat, Xf = Xf, Xm = Xm, type = type)
  W11.x <- W11.xtilde(xtilde = par, W11_mr1 = W11_mr1, Xf = Xf, Xm = Xm, uhat=  uhat, theta_m = theta_m, type = type)
  
  ## Dimension of nf + nm +1 x s
  dW11.dxtilde <- dW11(xtilde = xtilde, uhat = uhat, Xm = Xm, Xf = Xf, Xfuhat = Xfuhat, theta_m = theta_m, type = type)
  ## Dimension of nf x s
  dW12.dxtilde <- dW12(xtilde = xtilde, Xf = Xf, Xm = Xm, theta_b = theta_b, theta_m = theta_m, uhat = uhat, Xfuhat = Xfuhat, type = type)
  
  dW12.temp <- dW11.temp <- matrix(0, ncol = (nrow(Xf)+nrow(Xm)+1), nrow = (nrow(Xf)+ nrow(Xm)+1))
  
  dxtildeout <- rep(NA, times = length(xtilde))
  
  for(s in 1:length(xtilde)){
    dW11.temp[,ncol(dW11.temp)] <- dW11.dxtilde[s,]
    dW11.temp[nrow(dW11.temp), ] <- dW11.dxtilde[s,]
    dW12.temp[1:nrow(Xf),ncol(dW12.temp)] <- dW12.dxtilde[s,]
    
    dKi.s <- dKi.dxtilde[[s]]
    
    ## related to W11
    sumoftrace <- sum(dKi.s * (W11.x* Vm^2)) + sum(Kn1i * (dW11.temp* Vm^2))
    
    ## related to W12
    sumoftrace <- sumoftrace + sum(dKi.s *(W12.x * Vm*Vb)) + sum(dKi.s * (t(W12.x) * Vm*Vb)) + sum(Kn1i *(dW12.temp *Vm*Vb)) + sum(Kn1i * (t(dW12.temp) *Vm*Vb)) 
    ## related to W22
    
    sumoftrace <- sumoftrace + sum(dKi.s*(W22 *Vb*Vb)) 
    
    dxtildeout[s] <- -sumoftrace
  }
  
  if(log.IMSPE){
    temp <- sum(Kn1i * (W11.x * Vm^2)) + 2*sum(Kn1i *(W12.x*Vm * Vb)) + sum(Kn1i*(W22 * Vb^2))
    
    temp <- Vb + Vm*(1+eps) - temp
    
    dxtildeout <- dxtildeout/temp
    
    
  }
  
  return(dxtildeout)
  
}



dKi <- function(xtilde = NULL, uhat = NULL, Xm = NULL, Xf = NULL, Xfuhat = NULL, theta_m = NULL, Ki = NULL, Km = NULL, Vm = NULL, K_x.fm = NULL, b = NULL, Kikx = NULL,
                eps = NULL, type = NULL){
  
  dKi.dxtilde <- list()
  
  dKi.dxtilde.temp <- matrix(NA, ncol = ncol(Ki)+1, nrow = nrow(Ki) + 1)
  
  
  bi <- 1/b
  
  dA21 <- dA21.gaussian(xtilde = xtilde, Xf = Xf, Xfuhat = Xfuhat, Xm = Xm, uhat =  uhat, theta_m = theta_m)*Vm
  tKxKx <- t(K_x.fm) %*% K_x.fm
  KikxtKikx <- tcrossprod(Kikx)
  
  for(s in 1:length(xtilde)){
    
    da21.s <- dA21[s, , drop = FALSE]
    
    Kitda12.s <- solve(Km,t(da21.s))
    
    ## Derivative of Bi

    dbi <- bi^2*drop(da21.s %*% Kikx + t(Kikx) %*% t(da21.s))
    
    dKi.dxtilde.temp[nrow(Ki)+1, nrow(Ki)+1] <- dbi
    
    tl <-  (dbi* KikxtKikx + bi*(tcrossprod(Kitda12.s, Kikx) + tcrossprod(Kikx, Kitda12.s)))
    
    dKi.dxtilde.temp[1:nrow(Ki), 1:nrow(Ki)] <- tl
    
    tr <- -(dbi * Kikx + bi * Kitda12.s)
    dKi.dxtilde.temp[1:nrow(Ki),(nrow(Ki)+1)] <- tr
    dKi.dxtilde.temp[(nrow(Ki)+1), 1:nrow(Ki)] <- t(tr)
    
    
    dKi.dxtilde[[s]] <- dKi.dxtilde.temp
    
  }
  
  return(dKi.dxtilde)
  
}



dA21.gaussian <- function(xtilde = NULL, Xf = NULL, Xfuhat = NULL, Xm = NULL, uhat =  NULL, theta_m = NULL){
  nm <- nrow(Xm)
  nf <- nrow(Xf)
  
  A.mat <- dA.mat <- dA.out <-  matrix(NA, nrow = length(xtilde), ncol = nm + nf)
  
  XfXm <- rbind(Xf, Xm[,1:ncol(Xf), drop = FALSE])
  
  ## Matrix elements related to all data in X space
  
  for(i in 1:ncol(Xf)){
    X <- XfXm[,i]
    x <- xtilde[i]
    tm <- theta_m[i]
    
    a <- (X-x)
    da <- 2*a/tm
    a <- -(a^2)/tm
    da <- da*exp(a)
    
    A.mat[i,] <- a
    dA.mat[i, ] <- da
    
  }
  
  ## Elements of matricies related to field data in U space
  
  a <- (uhat - xtilde[(ncol(Xf)+1):ncol(Xm)])
  da <- 2*a/theta_m[(ncol(Xf)+1):ncol(Xm)]
  a <- -(a^2)/theta_m[(ncol(Xf)+1):ncol(Xm)]
  da <- da*exp(a)
  
  A.mat[(ncol(Xf)+1):ncol(Xm),1:nf] <- a
  dA.mat[(ncol(Xf)+1):ncol(Xm),1:nf] <- da
  
  ## Elements of matricies related to model data in U space
  
  for(i in (ncol(Xf)+1):ncol(Xm)){
    X <- Xm[,i]
    x <- xtilde[i]
    tm <- theta_m[i]
    
    a <- (X-x)
    da <- 2*a/tm
    a <- -(a^2)/tm
    da <- da*exp(a)
    
    A.mat[i,(nf +1):(nm + nf) ] <- a
    dA.mat[i,(nf+1):(nm+nf) ] <- da
    
  }
  
  for(i in 1:nrow(A.mat)){
    dA.out[i,] <- exp(colSums(A.mat[-i, , drop=  FALSE]))*dA.mat[i,]
  }
  
  return(dA.out)
  
}

dW11 <- function(xtilde = NULL, uhat = NULL, Xm = NULL, Xf = NULL, Xfuhat = NULL, theta_m = NULL, type = NULL){
  
  W.mat <- dW.mat <-  matrix(0, ncol = (nrow(Xf) + nrow(Xm) + 1), nrow = length(xtilde))
  
  Xfmb <- rbind(Xf, Xm[,1:ncol(Xf), drop = FALSE])
  
  
  ## Elements related to all data in X space
  W.mat[1:ncol(Xf),] <- t(apply(rbind(Xf,Xm[,1:ncol(Xf), drop = FALSE], xtilde[1:ncol(Xf)],theta_m[1:ncol(Xf)]), 2, function(X){
    x <- X[length(X)-1]
    tm <- X[length(X)]
    X <- X[-length(X)]
    
    w <- log(hetGP::Wij(matrix(X, ncol= 1, nrow = length(X)), mu2 = matrix(x, ncol = 1, nrow = 1), theta = tm, type = type))
    return(w)
    
  }))
  ## Apply below not necessary?
  
  ## Elements related to model data in U space - exponential difference between uhat and model locations
  W.mat[(ncol(Xf)+1):ncol(Xm), (nrow(Xf)+1):(nrow(Xm)+nrow(Xf))] <- t(apply(rbind(Xm[,(ncol(Xf)+1):ncol(Xm), drop = FALSE], uhat,theta_m[(ncol(Xf)+1):ncol(Xm)]), 2, function(X){
    u <- X[length(X)-1]
    tm <- X[length(X)]
    X <- X[-c(length(X)-1,length(X))]
    
    w <- -(X-u)^2/tm
    
    return(w)
  }))
  
  # exponentiated distance between uhat and xtilde
  uhatutilde <- -((xtilde[(ncol(Xf)+1):ncol(Xm)] - uhat)^2)/theta_m[(ncol(Xf)+1):ncol(Xm)]
  # Elements related to model data in U space
  W.mat[(ncol(Xf)+1):ncol(Xm), (nrow(Xf)+1):(nrow(Xm)+nrow(Xf))] <-   W.mat[(ncol(Xf)+1):ncol(Xm), (nrow(Xf)+1):(nrow(Xm)+nrow(Xf))] + uhatutilde
  
  # Elements related to field data in U space - kernel related to expontiated distance between xtilde and uhat
  W.mat[(ncol(Xf)+1):(ncol(Xm)), 1:nrow(Xf)] <- uhatutilde
  
  # Elements related to xtilde in U space
  W.mat[(ncol(Xf)+ 1):ncol(Xm), (nrow(Xm) + nrow(Xf) +1)] <- 2* uhatutilde
  
  # Elements related to model data in U space
  dW.mat[(ncol(Xf)+1):ncol(Xm), (nrow(Xf)+1):(nrow(Xm)+nrow(Xf))] <- 2*((uhat - xtilde[(ncol(Xf)+1):ncol(Xm)])/theta_m[(ncol(Xf)+1):ncol(Xm)])*exp(W.mat[(ncol(Xf)+1):ncol(Xm), (nrow(Xf)+1):(nrow(Xm)+nrow(Xf)), drop = FALSE])
  
  
  
  # Elements related to all locations in X space
  
  dW.mat[1:ncol(Xf),1:(nrow(Xm)+ nrow(Xf)+1)] <- t(apply(rbind(Xf,Xm[,1:ncol(Xf), drop = FALSE], xtilde[1:ncol(Xf)],theta_m[1:ncol(Xf)]), 2, function(X){
    x1<- X[length(X)-1]
    tm <- X[length(X)]
    X <- X[-c(length(X)-1, length(X))]
    
    w <- c11.hetGP.gaussian(x1 = x1, x2 = X, theta = tm)
    w <- c(w,exp(-2*(x1^2)/tm) - exp(-(2*(x1-1)^2)/(tm)))
    return(w)
    
  }))

  ## Elements related to field data in U space
  dW.mat[(ncol(Xf)+1):(ncol(Xm)), 1:nrow(Xf)] <-  2*((uhat - xtilde[(ncol(Xf)+1):ncol(Xm)])/theta_m[(ncol(Xf)+1):ncol(Xm)])*exp(uhatutilde)
  
  ## Elements related to xtilde in Uspace
  dW.mat[(ncol(Xf)+1):ncol(Xm), (nrow(Xf)+ nrow(Xm) + 1)] <- 4*((uhat - xtilde[(ncol(Xf)+1):ncol(Xm)])/theta_m[(ncol(Xf)+1):ncol(Xm)])*exp(2*uhatutilde)
  
  dW.out <- matrix(NA, ncol = nrow(Xm) + nrow(Xf) + 1, nrow = length(xtilde))
  
  
  for(i in 1:length(xtilde)){
    dW.out[i,] <- exp(colSums(W.mat[-i, , drop = FALSE]))*dW.mat[i,]
  }
  
  ## note, dW.out provides the row vector [c_1^T c_2] 
  return(dW.out)
  
}

dW12 <- function(xtilde = NULL, Xf = NULL, Xm = NULL, theta_b = NULL, theta_m = NULL, uhat = NULL, Xfuhat = NULL, type = NULL){
  
  dW.out <- W.mat <- dW.mat <- matrix(0, ncol = nrow(Xf), nrow = length(xtilde))
  
  W.mat[1:ncol(Xf),] <- t(apply(rbind(Xf, xtilde[1:ncol(Xf)], theta_m[1:ncol(Xf)], theta_b[1:ncol(Xf)]), 2, function(X){
    
    tb <- X[length(X)]
    tm <- X[length(X)-1]
    x <- X[length(X)-2]
    X <- X[-(length(X) - c(0,1,2))]
    
    a <- (X*tm + x*tb)/(tm + tb)
    b <- a -1
    sqrt.ls <- sqrt((1/tm + 1/tb)^(-1))
    
    w <- log(erf(a/sqrt.ls) - erf(b/sqrt.ls))
    
    w <- -(x - X)^2/(tm + tb) + log(0.5) + 0.5*log(pi) - 0.5*log(1/tm + 1/tb) + w
    
    return(w)
  }))
  
  W.mat[(ncol(Xf)+1):ncol(Xm),] <- -(uhat - xtilde[(ncol(Xf)+1):ncol(Xm)])^2/(theta_m[(ncol(Xf)+1):ncol(Xm)])
  
  ## elements related to X space
  
  dW.mat[1:ncol(Xf),] <- t(apply(rbind(Xf, xtilde[1:ncol(Xf)], theta_m[1:ncol(Xf)], theta_b[1:ncol(Xf)]),2, function(X){
    
    tb <- X[length(X)]
    tm <- X[length(X)-1]
    x <- X[length(X)-2]
    X <- X[-(length(X) - c(0,1,2))]
    
    w <- c22.gaussian(x1 = x, x2 = X, tm = tm, tb = tb)
    
    return(w)
  }))
  
  dW.mat[(ncol(Xf)+1):ncol(Xm),] <- 2*(uhat - xtilde[(ncol(Xf)+1):ncol(Xm)])/(theta_m[(ncol(Xf)+1):ncol(Xm)])*exp(-(uhat-xtilde[(ncol(Xf)+1):ncol(Xm)])^2/theta_m[(ncol(Xf)+1):ncol(Xm)])
 
  for(i in 1:length(xtilde)){
    dW.out[i,] <- exp(colSums(W.mat[-i, , drop = FALSE]))*dW.mat[i,]
  }
  
  return(dW.out)
}

c11.hetGP.gaussian <- function(x1 = NULL, x2 = NULL, theta = NULL){
  ## differentiation wrt x1
  
  arg1 <- sqrt(pi/2)*exp(-(x2-x1)^2/(2*theta))
  
  arg2 <- (x2-x1)*(erf((2-(x2+x1))/sqrt(2*theta)) + erf((x1+x2)/sqrt(2*theta)))/(2*sqrt(theta))
  
  arg3 <- 0.5*sqrt(2/pi)*(exp(-(x1+x2)^2/(2*theta)) - exp(-(2-x1-x2)^2/(2*theta)))
  
  dw <- arg1*(arg2 + arg3)
  
  return(dw)  
}

c22.gaussian <- function(x1 = NULL, x2 = NULL, tm = NULL, tb = NULL){
  
  ## differentiation wrt x1

  a <- (tm*x2 + tb*x1)/(tm + tb)
  b <- a - 1
  sqrt.ls <- sqrt((1/tm + 1/tb)^(-1))
  
  a <- a/sqrt.ls
  b <- b/sqrt.ls
  
  dw <- sqrt(pi)*sqrt.ls*(x2-x1)*(erf(a)-erf(b)) + tb*(exp(-(a^2)) - exp(-(b^2)))
  dw <- exp(-(x1-x2)^2/(tm + tb))/(tm + tb)*dw
  
  return(dw)
  
}

############ Miscelaneous #####################

## Error Function

erf <- function(x){2* pnorm(x * sqrt(2)) -1}

## Performance Functions

## yM is the response of the observed model data
## yF is the response of the observed field data
## Xf.test is the design of the testing set for the field data
## yF.true.test is the deterministic set of field data corresponding to Xf.test
## yF.noise.test is the noisy set of fiel data corresponding to Xf.test
## See the IMSPE.optim function in koh-imspe.R for additional variable definitions

rmse.score.fn <- function(Xm = NULL, yM = NULL, Xf = NULL, yF = NULL, uhat = NULL, Vm = NULL, Vb = NULL, theta_m = NULL, 
                          theta_b = NULL, gb = NULL, Xf.test = NULL, yF.true.test = NULL, eps = NULL, type = NULL){
  
  Xfm <- cbind(Xf, rep(1, times = nrow(Xf)) %x% t(uhat))
  Xf.test.uhat <- cbind(Xf.test, rep(1, times = nrow(Xf.test)) %x% t(uhat))
  
  Km <- cov_gen(Xm, theta = theta_m, type = type) + diag(nrow(Xm))*eps
  km.pred <- cov_gen(Xm, Xf.test.uhat, theta = theta_m, type = type)
  km.f <- cov_gen(Xm, Xfm, theta = theta_m, type = type)
  
  KmiyM <- solve(Km, yM)

  ## Predictive mean for testing set
  mu.m <- crossprod(km.pred, KmiyM)

  
  ## Predictive mean for field data given uhat
  
  Ym <- crossprod(km.f, KmiyM)
  
  ## Predicting GPs for field data seperatly for numerical stability
  
  Yb <- yF - Ym
  
  Kb <- cov_gen(Xf, theta = theta_b, type = type) + diag(nrow(Xf))*gb
  
  kb.pred <- cov_gen(Xf, Xf.test, theta = theta_b, type = type)
  
  KbiYb <- solve(Kb, Yb)
  
  mu.b <- crossprod(kb.pred, KbiYb)

  mu.pred <- mu.m + mu.b
  
  ### RMSE
  
  rmse <- sqrt(mean((yF.true.test -mu.pred)^2))

  return(list(rmse = rmse))
  
}

sx.score <- function(Xm = NULL, yM = NULL, Xf = NULL, yF = NULL, uhat = NULL, Vm = NULL, Vb = NULL, theta_m = NULL, 
                              theta_b = NULL, gb = NULL, Xf.test = NULL, yF.test = NULL, eps = NULL, type = NULL){
  
  Xuhatf <- cbind(Xf, rep(1, times = nrow(Xf)) %x% t(uhat))
  Xf.test.uhat <- cbind(Xf.test, rep(1, times = nrow(Xf.test)) %x% t(uhat))
  
  #### Full covariance predictive mean
  
  Xfm <- rbind(Xuhatf, Xm)
  
  Km.full <- cov_gen(Xfm, theta = theta_m, type = "Gaussian")
  diag(Km.full) <- diag(Km.full) + eps
  
  Kmb.full <- Km.full <- Vm*Km.full
  Kb <- cov_gen(Xf, theta = theta_b, type = "Gaussian")
  diag(Kb) <- diag(Kb) + gb
  Kmb.full[1:nrow(Xf), 1:nrow(Xf)] <-  Kmb.full[1:nrow(Xf), 1:nrow(Xf)] + Vb*Kb
  
  KmbiY <- solve(Kmb.full, c(yF, yM))
  
  kmbpred <- Vm*cov_gen(Xf.test.uhat, Xfm, theta = theta_m, type = "Gaussian")
  kmbpred[ ,1:nrow(Xf)] <- kmbpred[ ,1:nrow(Xf)] + Vb*cov_gen(Xf.test, Xf, theta = theta_b, type = "Gaussian")
  
  Yf.hat <- kmbpred %*% KmbiY
  
  mu.pred <- Yf.hat
  
  ### RMSE
  
  rmse <- sqrt(mean((yF.test -mu.pred)^2))
  
  ## Score
  
  nm <- nrow(Xm)
  nf <- nrow(Xf)
  
  k.pred <- Vm*cov_gen(X1 = Xfm, X2 = Xf.test.uhat, theta = theta_m, type = "Gaussian")
  k.pred[1:nf,] <- k.pred[1:nf,] + Vb*cov_gen(X1 = Xf, X2 = Xf.test, theta = theta_b, type = "Gaussian")
  
  kKik <- crossprod(k.pred,solve(Kmb.full, k.pred))
  
  Kmpred <- cov_gen(Xf.test.uhat, theta = theta_m, type = "Gaussian")
  diag(Kmpred) <- diag(Kmpred) + eps
  
  Kbpred <- cov_gen(Xf.test, theta = theta_b, type = "Gaussian")
  
  pred.var <- Vm*Kmpred + (Kbpred + diag(nrow(Kbpred))* gb)*Vb - kKik
  pred.var <- (pred.var + t(pred.var))/2
  
  if(any(diag(pred.var) < 0)){
    warning("Calculated negative predictive variance for score function", immediate. = TRUE)
  }
  
  bpred <- yF.test - mu.pred
  
  normal.kernel <- crossprod(bpred, solve(pred.var, bpred)) #(yF.test - mu.pred)^2/pred.var
  score <- -determinant(pred.var)$modulus - normal.kernel
  
  
  return(list(rmse = rmse, score = score))
}
