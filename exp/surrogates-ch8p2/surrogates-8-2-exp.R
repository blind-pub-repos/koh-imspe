#!/usr/bin/env Rscript

{
  if(Sys.getenv("RSTUDIO") == 1){
    setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
    expwd <- getwd()
  }else if(Sys.info()[4] == "mordor"){
    .libPaths(new = "~/Rmkl_libs")
    expwd <- getwd()
  }else{
    expwd <- getwd()
  }
  

source("surrogates-8-2-calib.R")
setwd("./../../")
source("code/koh-imspe.R")
setwd(expwd)


library(laGP)
library(plgp)
library(lhs)
library(hetGP)
library(geometry)
library(tgp)

library(profvis)

deleteGPseps()
deleteGPs()

#profvis({
eps <- sqrt(.Machine$double.eps)

mc.iters <- 100
min.pts <- 30
max.pts <- 130
tests <- c(min.pts, max.pts)
rndm.strts <- 600
testing.size <- 1000
field.sd <- 0.25
candidate <- FALSE
log.IMSPE <- TRUE
unique.field <- 25
field.reps <- 2
hyperparams.given <- FALSE
uhat.starts <- 10
sep <- TRUE

file.name <- "01-19-2023-surrogates-100reps.RData"

notes <- "fixed an error with the random design"

type <- "Gaussian"

test.types <- c("imspe", "lhs", "random")

param.rec <- matrix(NA, ncol=  mc.iters, nrow = max.pts)
 

performance.list <- list()
list.temp <- list()

list.temp$rmse <- param.rec
list.temp$uhat1 <- param.rec
list.temp$uhat2 <- param.rec

theta_m.list <- theta_b.list <- list()

for(i in 1:mc.iters){
  theta_m.list[[i]] <- matrix(NA, ncol = 4, nrow = max.pts)
  theta_b.list[[i]] <- matrix(NA, ncol = 2, nrow = max.pts)
}

list.temp$gb <- param.rec
list.temp$Vb <- param.rec
list.temp$Vm <- param.rec
list.temp$theta_m <- theta_m.list
list.temp$theta_b <- theta_b.list

for(i in test.types){
  performance.list[[i]] <- list.temp
}

## For each iteration of the monte carlo experiment

if(hyperparams.given){
  Vb <- 839.601060559771
  Vm <- 18.6286883891487
  theta_b1 <- 1.40053680891583 
  theta_b2 <- 1.47278070772745
  theta_m1 <- 0.0932946036561526
  theta_m2 <- 0.387151715458473 
  theta_m3 <- 6.17505379625791
  theta_m4 <- 5.74438473211422
  gb = 8.01598506407026e-05
  
  theta_b <- c(theta_b1, theta_b2)
  theta_m <- c(theta_m1, theta_m2, theta_m3, theta_m4)
}

for(i in 1:mc.iters){
  
  print(paste("mc iter", i, collapse = ""))
  
  ## Training Data
  

  
  Xm.full <- randomLHS(n= max.pts, k = 4)
  
  Xm <- Xm.full[1:min.pts, ,drop = FALSE]

  #Computer Model Evaluations
  yM <- Comp_Model(Xm[,1:2],Xm[,3:4])
  #Field Data
  Xf <- unname(as.matrix(expand.grid(seq(from  =0, to=  1, length.out = sqrt(unique.field)), seq(from= 0, to = 1, length.out = sqrt(unique.field)))))
  Xf <- rbind(Xf,Xf)
  
  #Field Observations
  yF <- drop(Field(Xf, sd = field.sd))
  
  ## Testing Data
  Xf.test <- randomLHS(n = testing.size, k = 2)
  yF.true.test <-  drop(Field(Xf.test,sd = 0))
  
  if(!hyperparams.given){
  previous.params <- NULL
  fits <- GP_Fits(Xm = Xm, yM = yM, Xf = Xf, yF = yF, eps = eps, starts = uhat.starts, sep = sep)
  
  initial.vals <- list(Xm = Xm, yM = yM, Xf = Xf, yF = yF, 
                       Vb = fits$tau2B, Vm = fits$tau2M, theta_m = fits$LS.M, theta_b = fits$LS.B, gb = fits$Nug.B, uhat = fits$uhat)
  
  }else{
    uhat <- uhat.est(Xm = Xm, Xf = Xf, yM = yM, yF = yF,eps = eps)
    initial.vals <- list(Xm = Xm, yM = yM, Xf = Xf, yF = yF, 
                         Vb = Vb, Vm = Vm, theta_m = theta_m, theta_b = theta_b, gb = gb, uhat = uhat)
  }
  
  rndm.tests <- matrix(runif(4*max.pts), ncol = ncol(Xm), nrow = max.pts)
  lhs.tests <- Xm.full[(min.pts +1):max.pts, , drop = FALSE]
  
  ## Iterate over scenarios
  for(s in test.types){
    
    
    print(paste(s, "design", collapse = ""))
    pb <- txtProgressBar(min =  min.pts, max = max.pts, initial = 0, style = 3)
    
    Xm <- initial.vals$Xm
    yM <- initial.vals$yM
    
    #Computer Model Evaluations
    
    Xf <- initial.vals$Xf
    yF <- initial.vals$yF
    
    Vb <- initial.vals$Vb
    Vm <- initial.vals$Vm
    theta_m <- initial.vals$theta_m
    theta_b <- initial.vals$theta_b
    gb <- initial.vals$gb
    uhat <- initial.vals$uhat
    
    previous.params <- list(theta_m = theta_m, theta_b = theta_b, gb = gb)
    
    current.perf <- rmse.score.fn(Xm = Xm, yM = yM, Xf = Xf, yF = yF, uhat = uhat, Vm = Vm, Vb = Vb, theta_m = theta_m,
                                  theta_b = theta_b, gb = gb, Xf.test = Xf.test, yF.true.test = yF.true.test,
                                  eps = eps)

    performance.list[[s]]$rmse[min.pts,i] <- current.perf$rmse
    performance.list[[s]]$uhat1[min.pts,i] <- uhat[1]
    performance.list[[s]]$uhat2[min.pts,i] <- uhat[2]
    performance.list[[s]]$theta_m[[i]][min.pts,] <- theta_m
    performance.list[[s]]$theta_b[[i]][min.pts,] <- theta_b
    performance.list[[s]]$gb[min.pts,i] <- gb
    performance.list[[s]]$Vb[min.pts,i] <- Vb
    performance.list[[s]]$Vm[min.pts,i] <- Vm

    
    ## Iterate over added points
    
    
    for(j in (min.pts+1):max.pts){
      
      if(s == "random"){
        xtilde <- rndm.tests[j,]
      }else if(s == "imspe"){
        xtilde <- IMSPE.optim(Xm = Xm, Xf = Xf, uhat = uhat, Vb = Vb, Vm = Vm,
                              theta_m = theta_m, theta_b = theta_b,gb = gb,
                              type = "Gaussian", starts = rndm.strts, eps = eps, candidate = candidate, 
                              log.IMSPE = log.IMSPE, optim.control = list(pgtol = 0.05, trace = 0, lmm = 13))
      }else if(s == "lhs"){
        xtilde <- lhs.tests[j-min.pts,]
      }
      
      Xm <- rbind(Xm,xtilde)
      
      yM <- c(yM, Comp_Model(Xm[j,1:2, drop = FALSE],Xm[j,3:4, drop = FALSE]))
      
      if(!hyperparams.given){
      fits <- GP_Fits(Xm = Xm, yM = yM, Xf = Xf, yF = yF, eps = eps, previous.params = previous.params, starts = uhat.starts, sep = sep)

      Vb <- fits$tau2B
      Vm <- fits$tau2M
      theta_m <- fits$LS.M
      theta_b <- fits$LS.B
      gb <- fits$Nug.B
      uhat <- fits$uhat
      
      previous.params <- list(theta_m = theta_m, theta_b = theta_b, gb = gb)
      }else{
        uhat <- uhat.est(Xm = Xm, Xf = Xf, yM = yM, yF = yF,eps = eps)
      }
      
      current.perf <- rmse.score.fn(Xm = Xm, yM = yM, Xf = Xf, yF = yF, uhat = uhat, Vm = Vm, Vb = Vb, theta_m = theta_m,
                                    theta_b = theta_b, gb = gb, Xf.test = Xf.test, yF.true.test = yF.true.test,
                                    eps = eps)
      
      performance.list[[s]]$rmse[j,i] <- current.perf$rmse
      performance.list[[s]]$uhat1[j,i] <- uhat[1]
      performance.list[[s]]$uhat2[j,i] <- uhat[2]
      performance.list[[s]]$theta_m[[i]][j,] <- theta_m
      performance.list[[s]]$theta_b[[i]][j,] <- theta_b
      performance.list[[s]]$gb[j,i] <- gb
      performance.list[[s]]$Vb[j,i] <- Vb
      performance.list[[s]]$Vm[j,i] <- Vm
      
      if(j == max.pts){
        performance.list[[s]]$data[[i]] <- list(Xm = Xm, Xf = Xf, yM = yM, yF = yF, yF.true.test = yF.true.test, Xf.test = Xf.test)
      }
      
      setTxtProgressBar(pb, value = j)
      
    }
    
    close(pb)
    
    
  }
  
  
  
  mcparams <- list(starts = rndm.strts, min.pts = min.pts, max.pts = max.pts,tests = tests, field.std = field.sd, field.dat = Xf, 
                   field.npoints = length(Xf), field.reps = 2, score.noise = TRUE, field.design.desc = "grid",  timestamp = Sys.time(),
                   notes = notes, hyperparams.given = hyperparams.given)
  

  performance.list$mcparams <- mcparams
  

  save(performance.list, file = file.name)
  
}


}