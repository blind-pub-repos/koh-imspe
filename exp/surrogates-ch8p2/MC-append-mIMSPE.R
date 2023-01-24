#!/usr/bin/env Rscript

######
### This R script adds a random design and simulator based IMSPE design to an
### existing MC experiment
#######
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
  
  library(PeriodicTable)
  library(laGP)
  library(plgp)
  library(lhs)
  library(hetGP)
  library(geometry)
  library(tgp)
  
  library(profvis)
  
  eps <- sqrt(.Machine$double.eps)
  
  ######### rework of hetGP:::IMSPE_optim
  
 
IMSPE.search <- function (model, replicate = FALSE, Xcand = NULL, control = list(tol_dist = 1e-06, 
                                                                                   tol_diff = 1e-06, multi.start = 20, maxit = 100, maximin = TRUE, 
                                                                                   Xstart = NULL), Wijs = NULL, seed = NULL, ncores = 1) 
  {
    if (replicate) {
      res <- unlist(mclapply(1:nrow(model$X0), crit_IMSPE, 
                             Wijs = Wijs, model = model, x = NULL, mc.cores = ncores))
      return(list(par = model$X0[which.min(res), , drop = FALSE], 
                  value = min(res), new = FALSE, id = which.min(res)))
    }
    if (is.null(control)) 
      control <- list(multi.start = 20, maxit = 100)
    if (is.null(control$multi.start)) 
      control$multi.start <- 20
    if (is.null(control$maxit)) 
      control$maxit <- 100
    if (is.null(control$maximin)) 
      control$maximin <- TRUE
    if (is.null(control$tol_dist)) 
      control$tol_dist <- 1e-06
    if (is.null(control$tol_diff)) 
      control$tol_diff <- 1e-06
    d <- ncol(model$X0)
    if (is.null(Wijs)) 
      Wijs <- Wij(mu1 = model$X0, theta = model$theta, type = model$covtype)
    if (is.null(Xcand)) {
      if (!is.null(control$Xstart)) {
        Xstart <- control$Xstart
      }
      else {
        if (is.null(seed)) 
          seed <- sample(1:2^15, 1)
        if (control$maximin) {
          if (d == 1) {
            Xstart <- matrix(seq(1/2, control$multi.start - 
                                   1/2, length.out = control$multi.start) + 
                               runif(control$multi.start, min = -1/2, max = 1/2), 
                             ncol = 1)/control$multi.start
          }
          else {
            Xstart <- maximinSA_LHS(lhsDesign(control$multi.start, 
                                              d, seed = seed)$design)$design
          }
        }
        else {
          Xstart <- lhsDesign(control$multi.start, d, seed = seed)$design
        }
      }
      res <- list(par = NA, value = Inf, new = NA)
      local_opt_fun <- function(i) {
        out <- try(optim(Xstart[i, , drop = FALSE], crit_IMSPE, 
                         method = "L-BFGS-B", lower = rep(0, d), upper = rep(1, 
                                                                             d), Wijs = Wijs, model = model, control = list(maxit = control$maxit), 
                         gr = deriv_crit_IMSPE))
        if (class(out) == "try-error") 
          return(NULL)
        
        if(out$convergence %in% c(51,52)){
          warning(paste("Optim for IMSPE failed to converge:  Code", out$convergence, sep = " ", collapse = " "), immediate. = TRUE)
        }
        return(out)
      }
      all_res <- mclapply(1:nrow(Xstart), local_opt_fun, mc.cores = ncores)
      res_min <- which.min(Reduce(c, lapply(all_res, function(x) x$value)))
      res <- list(par = apply(all_res[[res_min]]$par, c(1, 
                                                        2), function(x) max(min(x, 1), 0)), value = all_res[[res_min]]$value, 
                  new = TRUE, id = NULL)
      if (control$tol_dist > 0 && control$tol_diff > 0) {
        dists <- sqrt(distance_cpp(res$par, model$X0))
        if (min(dists) < control$tol_dist) {
          res <- list(par = model$X0[which.min(dists), 
                                     , drop = FALSE], value = crit_IMSPE(x = model$X0[which.min(dists), 
                                                                                      , drop = F], model = model, id = which.min(dists), 
                                                                         Wijs = Wijs), new = FALSE, id = which.min(dists))
        }
        else {
          id_closest <- which.min(dists)
          imspe_rep <- crit_IMSPE(model = model, id = id_closest, 
                                  Wijs = Wijs)
          if ((imspe_rep - res$value)/res$value < control$tol_diff) {
            res <- list(par = model$X0[which.min(dists), 
                                       , drop = FALSE], value = imspe_rep, new = FALSE, 
                        id = which.min(dists))
          }
        }
      }
      return(res)
    }
    else {
      crit_IMSPE_mcl <- function(i, model, Wijs, Xcand) {
        crit_IMSPE(x = Xcand[i, , drop = F], model = model, 
                   Wijs = Wijs)
      }
      res <- unlist(mclapply(1:nrow(Xcand), crit_IMSPE_mcl, 
                             Xcand = Xcand, Wijs = Wijs, model = model, mc.cores = ncores))
      tmp <- which(duplicated(rbind(model$X0, Xcand[which.min(res), 
                                                    , drop = FALSE]), fromLast = TRUE))
      if (length(tmp) > 0) 
        return(list(par = Xcand[which.min(res), , drop = FALSE], 
                    value = min(res), new = FALSE, id = tmp))
      return(list(par = Xcand[which.min(res), , drop = FALSE], 
                  value = min(res), new = TRUE, id = NULL))
    }
  }
  
source("surrogates-8-2-calib.R")
setwd("./../../")
source("code/koh-imspe.R")
setwd(expwd)
  
  previous.file <- "04-07-2022-surrogates-100reps.RData"
  previous.file.short <- "04-07-2022-surrogates-100reps"
  file.name <- previous.file
  
  load(previous.file)
  previous.data <- performance.list
  
  mc.iters <- 100
  min.pts <- previous.data$mcparams$min.pts
  max.pts <- previous.data$mcparams$max.pts
  tests <- c(min.pts, max.pts)
  imspe.strts <- previous.data$mcparams$starts
  candidate <- FALSE
  log.IMSPE <- TRUE
  hyperparams.given <- FALSE

  deleteGPseps()
  deleteGPs()
  
  p <- 2
  d <- 4
  
  type <- "Gaussian"
  
  file.name <- paste(previous.file.short,"-appended.RData", collapse = "", sep = "")
  
  notes <- "attempting to append additional tests"
  
  test.types <- c("m-imspe")
  
  param.rec <- matrix(NA, ncol=  mc.iters, nrow = max.pts)
  
  
  #performance.list <- list()
  list.temp <- list()
  
  list.temp$rmse <- param.rec
  list.temp$uhat1 <- param.rec
  list.temp$uhat2 <- param.rec
  
  theta_m.list <- theta_b.list <- list()
  
  for(i in 1:mc.iters){
    theta_m.list[[i]] <- matrix(NA, ncol = d, nrow = max.pts)
    theta_b.list[[i]] <- matrix(NA, ncol = p, nrow = max.pts)
  }
  
  list.temp$gb <- param.rec
  list.temp$Vb <- param.rec
  list.temp$Vm <- param.rec
  list.temp$theta_m <- theta_m.list
  list.temp$theta_b <- theta_b.list
  
  for(i in test.types){
    performance.list[[i]] <- list.temp
  }
  
  for(i in 1:mc.iters){
    
    #Field Data
    Xf <- previous.data$imspe$data[[i]]$Xf
    Xf.test <- previous.data$imspe$data[[i]]$Xf.test
    
    #Field Observations
    
    yF <- previous.data$imspe$data[[i]]$yF
    yF.true.test <- previous.data$imspe$data[[i]]$yF.true.test
    
    print(paste("MC Run ", i, collapse = ""))
    
    Xm <- previous.data$imspe$data[[i]]$Xm[1:min.pts, , drop = FALSE]
    yM <- previous.data$imspe$data[[i]]$yM[1:min.pts]
    
    
    if("random" %in% test.types){
      rndm.tests <- matrix(runif(d*max.pts), ncol = ncol(Xm), nrow = max.pts)
      yM.rndm <- rep(NA, times= max.pts)
      print("Gathering computer sim data for random design in parallel")
      yM.rndm[(min.pts+1):max.pts] <- Comp_Model(rndm.tests[(min.pts + 1):max.pts,1:p],rndm.tests[(min.pts + 1):max.pts,(p + 1):d], ranges = ranges, t = t, h = h, outs = 2)
      print("Done!")
    }
    
    
    previous.params <- NULL
    print("Fitting initial GP")
    fits <- list(gb = previous.data$imspe$gb[min.pts,i], Vb = previous.data$imspe$Vb[min.pts,i], Vm = previous.data$imspe$Vm[min.pts,i],
                 theta_m = previous.data$imspe$theta_m[[i]][min.pts,], theta_b = previous.data$imspe$theta_b[[i]][min.pts,],
                 uhat = c(previous.data$imspe$uhat1[min.pts,i], previous.data$imspe$uhat2[min.pts,i]))

    print("Done!")
    
    initial.vals <- list(Xm = Xm, yM = yM, Xf = Xf, yF = yF, 
                         Vb = fits$Vb, Vm = fits$Vm, theta_m = fits$theta_m, theta_b = fits$theta_b, gb = fits$gb, uhat = fits$uhat)
    
    
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
      performance.list[[s]]$score[min.pts,i] <- current.perf$score
      performance.list[[s]]$uhat1[min.pts,i] <- uhat[1]
      performance.list[[s]]$uhat2[min.pts,i] <- uhat[2]
      performance.list[[s]]$theta_m[[i]][min.pts,] <- theta_m
      performance.list[[s]]$theta_b[[i]][min.pts,] <- theta_b
      performance.list[[s]]$gb[min.pts,i] <- gb
      performance.list[[s]]$Vb[min.pts,i] <- Vb
      performance.list[[s]]$Vm[min.pts,i] <- Vm
      
      ## Iterate over added points
      
      for(j in (min.pts+1):max.pts){
        print(paste("finding", s,  "point", collapse = " "))
        if(s == "random"){
          xtilde <- rndm.tests[j,]
          Xm <- rbind(Xm,xtilde)
          
          yM <- c(yM, yM.rndm[j])
        }else if(s == "m-imspe"){
          
          hgp <- mleHomGP(X = Xm, Z = yM, known = list(theta = theta_m, g = eps), covtype = "Gaussian",
                          maxit = 0, settings = list(return.Ki = TRUE))
          
          xtilde <- IMSPE_optim(model = hgp, h = 0)$par
          set.seed(NULL)
          
          Xm <- rbind(Xm,xtilde)
          
          yM <- c(yM, Comp_Model(Xm[j,1:p, drop = FALSE],Xm[j,(p +1):d, drop = FALSE]))
        }
        
        print("fitting new GP")
        fits <- GP_Fits(Xm = Xm, yM = yM, Xf = Xf, yF = yF, eps = eps, previous.params = previous.params)
        
        
        Vb <- fits$tau2B
        Vm <- fits$tau2M
        theta_m <- fits$LS.M
        theta_b <- fits$LS.B
        gb <- fits$Nug.B
        uhat <- fits$uhat
        
        previous.params <- list(theta_m = theta_m, theta_b = theta_b, gb = gb)
        
        current.perf <- rmse.score.fn(Xm = Xm, yM = yM, Xf = Xf, yF = yF, uhat = uhat, Vm = Vm, Vb = Vb, theta_m = theta_m,
                                      theta_b = theta_b, gb = gb, Xf.test = Xf.test, yF.true.test = yF.true.test,
                                      eps = eps)
        
        
        performance.list[[s]]$rmse[j,i] <- current.perf$rmse
        performance.list[[s]]$score[j,i] <- current.perf$score
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
    
    test.types.plot <- c("random", "lhs", "imspe", "m-imspe")
    
    pts.important <- min.pts:max.pts
    run.vec <- 1:i
    
    imspe.plt <- apply(performance.list$imspe$rmse[pts.important,run.vec, drop = FALSE],1,mean)
    plt.ylim <- range(imspe.plt)
    if("random" %in% test.types.plot){
      random.plt <- apply(performance.list$random$rmse[pts.important,run.vec, drop = FALSE],1,mean)
      plt.ylim <- range(c(plt.ylim, random.plt))
    }
    if("lhs" %in% test.types.plot){
      lhs.plt <- apply(performance.list$lhs$rmse[pts.important,run.vec, drop = FALSE],1,mean)
      plt.ylim <- range(c(plt.ylim, lhs.plt))
    }
    if("m-imspe" %in% test.types.plot){
      m.imspe.plt <- apply(performance.list$`m-imspe`$rmse[pts.important,run.vec, drop = FALSE],1,mean)
      plt.ylim <- range(c(plt.ylim, m.imspe.plt))
    }
    
    plot(imspe.plt, type = "l", 
         col = "blue", main = paste("average behavior after run no", i, sep = " ", collapse = " "), ylab = "RMSE", ylim = plt.ylim)
    if("random" %in% test.types.plot){
      lines(random.plt, col = "red")
    }
    if("lhs" %in% test.types.plot){
      lines(lhs.plt, col = "black")
    }
    if("m-imspe" %in% test.types.plot){
      lines(m.imspe.plt, col = "green")
    }
    
    
    run.vec <- i
    
    imspe.plt <- apply(performance.list$imspe$rmse[pts.important,run.vec, drop = FALSE],1,mean)
    plt.ylim <- range(imspe.plt)
    if("random" %in% test.types.plot){
      random.plt <- apply(performance.list$random$rmse[pts.important,run.vec, drop = FALSE],1,mean)
      plt.ylim <- range(c(plt.ylim, random.plt))
    }
    if("lhs" %in% test.types.plot){
      lhs.plt <- apply(performance.list$lhs$rmse[pts.important,run.vec, drop = FALSE],1,mean)
      plt.ylim <- range(c(plt.ylim, lhs.plt))
    }
    if("m-imspe" %in% test.types.plot){
      m.imspe.plt <- apply(performance.list$`m-imspe`$rmse[pts.important,run.vec, drop = FALSE],1,mean)
      plt.ylim <- range(c(plt.ylim, m.imspe.plt))
    }
    
    plot(imspe.plt, type = "l", 
         col = "blue", main = paste("behavior of run no", i, sep = " ", collapse = " "), ylab = "RMSE", ylim = plt.ylim)
    if("random" %in% test.types.plot){
      lines(random.plt, col = "red")
    }
    if("lhs" %in% test.types.plot){
      lines(lhs.plt, col = "black")
    }
    if("m-imspe" %in% test.types.plot){
      lines(m.imspe.plt, col = "green")
    }

    
    
    mcparams <- list(imspe.starts = imspe.strts, min.pts = min.pts, max.pts = max.pts,tests = tests, score.noise = TRUE,  timestamp = Sys.time(),
                     notes = notes, hyperparams.given = hyperparams.given, mc.iters = mc.iters, test.types = test.types,
                     candidate = candidate)
    
    
    performance.list$mcparams <- mcparams
    
    
    
    
    save(performance.list, file = file.name)
    
    
    
  }
  
}