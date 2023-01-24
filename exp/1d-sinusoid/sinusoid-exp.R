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
  

source("sinusoid-calib.R")
setwd("./../../")
source("code/koh-imspe.R")
setwd(expwd)


library(laGP)
library(plgp)
library(lhs)
library(hetGP)
library(geometry)
library(tgp)

eps <- sqrt(.Machine$double.eps)

mc.iters <- 1000
min.pts <- 10
max.pts <- 50
tests <- c(min.pts,max.pts)
rndm.strts <- 500
field.sd <- 0.1
candidate <- FALSE
log.IMSPE <- TRUE
unique.field <- 5
field.reps <- 2

deleteGPs()
deleteGPseps()

## name for the output data file

file.name <- "11-18-2022-sinusoid-1000reps.RData"

notes <- "code submission version"

type <- "Gaussian"

############################
#####  in the "test.types" character vector "imspe" = KOH-IMSPE, "lhs" = LHS, 
#####  "random" = independent random uniform points, and "m-imspe" = M-IMSPE
#############################

test.types <- c("imspe", "lhs", "random", "m-imspe")

param.rec <- matrix(NA, ncol=  mc.iters, nrow = max.pts)

performance.list <- list()
list.temp <- list()

list.temp$rmse <- param.rec
list.temp$uhat <- param.rec

theta_m.list <- theta_b.list <- list()

for(i in 1:mc.iters){
  theta_m.list[[i]] <- matrix(NA, ncol = 2, nrow = max.pts)
  theta_b.list[[i]] <- matrix(NA, ncol = 1, nrow = max.pts)
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



for(i in 1:mc.iters){
  
  print(paste("mc iter", i, collapse = ""))
  
  ## Training Data
  
  model.pts <- min.pts
 
  Xm.full <- randomLHS(n= max.pts, k = 2)
  
  Xm <- Xm.full[1:model.pts, ,drop = FALSE]
  
  #Computer Model Evaluations
  yM <- Comp_Model(Xm[,1],Xm[,2])
  #Field Data
  Xf <- matrix(seq(from = 0, to = 1, length.out = unique.field), ncol= 1, nrow=  unique.field)
  Xf <- rbind(Xf,Xf)
  
  #Field Observations
  yF <- drop(Field(Xf, sd = field.sd))
  
  ## Testing Data
  Xf.test <- randomLHS(n = 100, k = 1)
  yF.true.test <-  drop(Field(Xf.test,sd = 0))
  

  previous.params <- NULL
  fits <- GP_Fits(Xm = Xm, yM = yM, Xf = Xf, yF = yF, eps = eps)
  initial.vals <- list(Xm = Xm, yM = yM, Xf = Xf, yF = yF,Vb = fits$tau2B, 
                       Vm = fits$tau2M, theta_m = fits$LS.M, theta_b = fits$LS.B, 
                       gb = fits$Nug.B, uhat = fits$uhat)

  
  
  rndm.tests <- matrix(runif(2*max.pts), ncol = 2, nrow = max.pts)
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
    performance.list[[s]]$uhat[min.pts,i] <- uhat
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
                              type = "Gaussian", p = 1, starts = rndm.strts, eps = eps, 
                              candidate = candidate, log.IMSPE = log.IMSPE,
                              optim.control = list(pgtol = "adaptive"))
      }else if(s == "lhs"){
        xtilde <- lhs.tests[j-min.pts,]
      }else if(s == "m-imspe"){
        hgp <- mleHomGP(X = Xm, Z = yM, known = list(theta = theta_m, g = eps), covtype = "Gaussian",
                        maxit = 0, settings = list(return.Ki = TRUE))
        
        xtilde <- IMSPE_optim(model = hgp, h = 0)$par
        set.seed(NULL)
      }
      
      Xm <- rbind(Xm,xtilde)
      rownames(Xm) <- NULL
      
      yM <- Comp_Model(Xm[,1],Xm[,2])
      yM <- unname(yM)
      
      fits <- GP_Fits(Xm = Xm, yM = yM, Xf = Xf, yF = yF, eps = eps, previous.params = previous.params)
      
      
      Vb <- fits$tau2B
      Vm <- fits$tau2M
      theta_m <- fits$LS.M
      theta_b <- fits$LS.B
      gb <- fits$Nug.B
      uhat <- fits$uhat
      
      previous.params <- list(theta_m = theta_m, theta_b = theta_b, gb = gb)
    
      
      current.perf <- rmse.score.fn(Xm = Xm, yM = yM, Xf = Xf, yF = yF, uhat = uhat, 
                                    Vm = Vm, Vb = Vb, theta_m = theta_m, theta_b = theta_b, 
                                    gb = gb, Xf.test = Xf.test, yF.true.test = yF.true.test,
                                    eps = eps)
      
      performance.list[[s]]$rmse[j,i] <- current.perf$rmse
      performance.list[[s]]$uhat[j,i] <- uhat
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
                   notes = notes)
  performance.list$mcparams <- mcparams
  #
  
  #
  save(performance.list, file = file.name)
  
  
  
}


}

