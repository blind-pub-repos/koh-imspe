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
  
  #############################################################################
  ### BEGIN Field data organization
  #############################################################################
  
  source("load-o-to-a-equal-1.R")
  
  
  dat <- read.csv("OA-0-to-one.csv")
  dat1 <- dat[,1:(which(names(dat)== "Sample.ID")-1)]
  dat2 <- dat[,(which(names(dat)== "Sample.ID")):(length(names(dat)))]
  
  m <- match(dat2$Sample.ID, dat1$Sample)
  dat1 <- dat1[m, ]
  
  
  dat <- cbind(dat1, dat2)
  
  dat$test <- NA
  dat$test[dat$Sample.Type == "Aqueous"] <- 1:(length(dat$test[dat$Sample.Type == "Aqueous"]))
  
  feeds <- which(dat$Sample.Type == "Feed")
  dat.list <- list()
  
  for(j in feeds){
    if(j == feeds[length(feeds)]){
      dat.temp <- dat[j:nrow(dat),]
    }else{
      dat.temp <- dat[j:(feeds[which(feeds == j) + 1]-1),]
    }
    for(i in dat.temp$test[!is.na(dat.temp$test)]){
      if(i == min(dat.temp$test, na.rm = TRUE)){
        tmp.var <- i
        next
        
      }else{
        dat.temp <- rbind(dat.temp[1:(which(dat.temp$test == i) - 1),], dat.temp[1,], dat.temp[(which(dat.temp$test == i)):nrow(dat.temp),])
        tmp.var <- c(tmp.var, i)
      }
      
    }  
    dat.temp$test <- unique(dat.temp$test[!is.na(dat.temp$test)]) %x% c(1,1)
    
    dat.list[[which(feeds == j)]] <- dat.temp
    
    
  }
  
  dat <- do.call("rbind", dat.list)
  dat$init.pH <- NA
  
  dat$init.pH[dat$Sample.Type == "Feed"] <- 1.99
  
  
  rmv.unit <- names(dat)[which(names(dat)== "Org.In"):which(names(dat) == "NaOH.Used")]
  
  for(i in rmv.unit){
    splt <- strsplit(dat[[i]],split = c("\\s+", ""))
    
    splt <- lapply(splt, function(X){
      X <- X[!(X == " " | X == "g" | X == "m" | X == "l")]
      X <- paste(X, sep = "", collapse = "")
      X <- as.numeric(X)
      return(X)
    })
    
    splt <- unlist(splt)
    
    dat[[i]] <- splt
    
    
  }
  
  ### Delete First test set due to inconsistent recording of NaOH usage
  
  dat <- dat[-(1:4),]
  
  dat[dat == "ND"] <- eps
  dat[dat == ""] <- NA
  
  ## Copnverting ug/L into mol/L
  
  elements <- c("Li", "Be", "Na", "Mg", "Al", "Si", "P", 
                "K", "Ca", "Sc", "Ti", "Cr", "Fe", "Mn", "Co", "Ni", "Cu", "Zn", 
                "Ga", "Rb", "Sr", "Y", "Nb", "Ba", "La", "Ce", "Pr", "Nd", "Sm", 
                "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Th", "U")
  
  for(i in elements){
    if(is.numeric(dat[[i]])){
      dat[[i]] <- dat[[i]]/1000000*mass(i)
    }else{
      splt <- strsplit(dat[[i]], split = "")
      splt <- lapply(splt, function(X){
        X <- X[!(X == ",")]
        X <- paste(X, sep = "", collapse = "")
        return(X)
      })
      splt <- unlist(splt)
      dat[[i]] <- (as.numeric(splt)/1000000)*mass(i)
    }
  }
  
  ### Add column for vol and mol of NaOH used
  ### Assume 
  ####  density of 50% NaOH is 1.52 g/cm^3
  ####  molarity of 50% NaOH is 18.9375 mol/L
  ####  density of 5% NaOH is 1.0554 g/cm^3
  ####  molarity of 5% NaOH is 1.32 mol/L
  
  
  dat$NaOH.vol <- NA
  dat$NaOH.mol <- NA
  
  dat$NaOH.vol[dat$Sample.Type == "Feed"] <- (dat$NaOH.Used[dat$Sample.Type == "Aqueous"]/1.52)/dat$Aq.In[dat$Sample.Type == "Aqueous"]
  
  dat$NaOH.mol[dat$Sample.Type == "Feed"] <- (dat$NaOH.Used[dat$Sample.Type == "Aqueous"]/1.52*18.9375)/dat$Aq.In[dat$Sample.Type == "Aqueous"]
  
  dat$OA <- NA
  dat$OA[which(!(is.na(dat$Org.In)))-1] <- (dat$Org.In/dat$Aq.In)[which(!(is.na(dat$Org.In)))]
  dat$OA <- log10(dat$OA)
  
  dat <- dat[c("test", "Na", "La","NaOH.vol", "NaOH.mol", "Sample.Type", "OA")]
  
  X <- dat[dat$Sample.Type == "Feed",]
  Y <- dat[!(dat$Sample.Type == "Feed"),]
  
  X$Sample.Type <- NULL
  X$test <- NULL
  
  Y$test <- NULL
  Y$NaOH.vol <- NULL
  Y$NaOH.mol <- NULL
  Y$Sample.Type <- NULL
  Y$OA <- NULL
  
  X <- rbind(X, X.1to1)
  Y <- rbind(Y, Y.1to1)
  
  mean.init.conc <- apply(X[1:2], 2, mean)
  mean.init.conc <- as.data.frame(matrix(mean.init.conc, ncol = length(mean.init.conc), nrow = 1))
  names(mean.init.conc) <- c("Na", "La")
  
  X$Na <- NULL
  X$La <- NULL
  
  ranges <- as.data.frame(apply(X, 2,range))
  #Structure of U:  [Na forward reaction, Na backward reaction, La forward reaction, La backward reaction]
  ranges <- cbind(ranges, data.frame(Na.f = c(-3,3), Na.b = c(-3,3), La.f = c(-3,3), La.b = c(-3,3)))
  ranges$OA <- c(-1,0)
  
  # ranges$Na[2] <- ranges$Na[2] + 0.1
  # ranges$Na[1] <- 0
  # ranges$La[2] <- ranges$La[2] + 0.1
  # ranges$La[1] <- 0
  
  X.code <- X
  
  for(j in 1:ncol(X)){
    X.code[,j] <- (X.code[,j] - ranges[1,j])/(ranges[2,j]-ranges[1,j])
  }
  
  Y.log <- log(Y)
  
  ###############################################################################################################
  ###############################################################################################################
  ###############################################################################################################
  ###############################################################################################################
  #############################################################################
  ### END Field data organization
  #############################################################################
  ###############################################################################################################
  ###############################################################################################################
  ###############################################################################################################
  ###############################################################################################################
  
  
  
  source("sx-sim copy.R")
  source("sx-calib.R")
  setwd("./../../")
  source("code/koh-imspe.R")
  setwd(expwd)
  #source("code/koh-imspe-nobias.R")
  
  previous.file <- "2022-11-07-sx-200-reps.RData"
  previous.file.short <- "2022-11-07-sx-200-reps"
  file.name <- previous.file
  
  load(previous.file)
  previous.data <- performance.list
  
  mc.iters <- previous.data$mcparams$mc.iters
  min.pts <- previous.data$mcparams$min.pts
  max.pts <- previous.data$mcparams$max.pts
  tests <- c(min.pts, max.pts)
  imspe.strts <- previous.data$mcparams$starts
  candidate <- FALSE
  log.IMSPE <- TRUE
  hyperparams.given <- FALSE
  ## uhat.starts <- 10 for the run starting 11/07/2022
  uhat.starts <- 10
  ## sep <- TRUE for the run starting on 11/07/2022
  sep <- TRUE
  deleteGPseps()
  deleteGPs()
  
  formals(Comp_Model)$log <- TRUE
  formals(RE.conc)$parallel <- TRUE
  formals(RE.conc)$init.conc <- mean.init.conc
  
  t <- previous.data$mcparams$RK.vars$t 
  h <- previous.data$mcparams$RK.vars$h
  
  p <- 3
  d <- 7
  
  type <- "Gaussian"
  
  file.name <- paste(previous.file.short,"-appended.RData", collapse = "", sep = "")
  
  notes <- "attempting to append additional tests"
  
  test.types <- c("random", "m-imspe")
  
  param.rec <- matrix(NA, ncol=  mc.iters, nrow = max.pts)
  
  
  #performance.list <- list()
  list.temp <- list()
  
  list.temp$rmse <- param.rec
  list.temp$score <- param.rec
  list.temp$uhat1 <- param.rec
  list.temp$uhat2 <- param.rec
  list.temp$uhat3 <- param.rec
  list.temp$uhat4 <- param.rec
  
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
    
    #hold.out <- sample(1:nrow(Y.log), size = floor(0.2*nrow(Y.log)))
    
    #Field Data
    Xf <- previous.data$imspe$data[[i]]$Xf
    Xf.test <- previous.data$imspe$data[[i]]$Xf.test
    #Xf <- unname(as.matrix(X.code[-hold.out,]))
    #Xf.test <- unname(as.matrix(X.code[hold.out, ,drop = FALSE]))
    
    #Field Observations
    
    yF <- previous.data$imspe$data[[i]]$yF
    yF.test <- previous.data$imspe$data[[i]]$yF.test
    
    #yF <- unname(Y.log[-hold.out,2,drop = TRUE])
    #yF.test <- unname(Y.log[hold.out,2,drop = TRUE])  #unname(Y.log[i,2,drop = TRUE])
    
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
                 uhat = c(previous.data$imspe$uhat1[min.pts,i], previous.data$imspe$uhat2[min.pts,i],
                          previous.data$imspe$uhat3[min.pts,i], previous.data$imspe$uhat4[min.pts,i]))
    #fits <- GP_Fits(Xm = Xm, yM = yM, Xf = Xf, yF = yF, eps = eps, starts = uhat.starts, sep = sep)
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
      
      current.perf <- sx.score(Xm = Xm, yM = yM, Xf = Xf, yF = yF, uhat = uhat, Vm = Vm, Vb = Vb, theta_m = theta_m,
                               theta_b = theta_b, gb = gb, Xf.test = Xf.test, yF.test = yF.test,
                               eps = eps)
      
      performance.list[[s]]$rmse[min.pts,i] <- current.perf$rmse
      performance.list[[s]]$score[min.pts,i] <- current.perf$score
      performance.list[[s]]$uhat1[min.pts,i] <- uhat[1]
      performance.list[[s]]$uhat2[min.pts,i] <- uhat[2]
      performance.list[[s]]$uhat3[min.pts,i] <- uhat[3]
      performance.list[[s]]$uhat4[min.pts,i] <- uhat[4]
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
          
          yM <- c(yM, Comp_Model(Xm[j,1:p, drop = FALSE],Xm[j,(p +1):d, drop = FALSE], ranges = ranges, t = t, h = h, outs = 2))
        }
        
        print("fitting new GP")
        fits <- GP_Fits(Xm = Xm, yM = yM, Xf = Xf, yF = yF, eps = eps, previous.params = previous.params, starts = uhat.starts, sep = sep)
        
        
        Vb <- fits$tau2B
        Vm <- fits$tau2M
        theta_m <- fits$LS.M
        theta_b <- fits$LS.B
        gb <- fits$Nug.B
        uhat <- fits$uhat
        
        previous.params <- list(theta_m = theta_m, theta_b = theta_b, gb = gb)
        
        current.perf <- sx.score(Xm = Xm, yM = yM, Xf = Xf, yF = yF, uhat = uhat, Vm = Vm, Vb = Vb, theta_m = theta_m,
                                 theta_b = theta_b, gb = gb, Xf.test = Xf.test, yF.test = yF.test,
                                 eps = eps)
        
        
        performance.list[[s]]$rmse[j,i] <- current.perf$rmse
        performance.list[[s]]$score[j,i] <- current.perf$score
        performance.list[[s]]$uhat1[j,i] <- uhat[1]
        performance.list[[s]]$uhat2[j,i] <- uhat[2]
        performance.list[[s]]$uhat3[j,i] <- uhat[3]
        performance.list[[s]]$uhat4[j,i] <- uhat[4]
        performance.list[[s]]$theta_m[[i]][j,] <- theta_m
        performance.list[[s]]$theta_b[[i]][j,] <- theta_b
        performance.list[[s]]$gb[j,i] <- gb
        performance.list[[s]]$Vb[j,i] <- Vb
        performance.list[[s]]$Vm[j,i] <- Vm
        
        if(j == max.pts){
          performance.list[[s]]$data[[i]] <- list(Xm = Xm, Xf = Xf, yM = yM, yF = yF, yF.test = yF.test, Xf.test = Xf.test)
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
    
    imspe.plt <- apply(performance.list$imspe$score[pts.important,run.vec, drop = FALSE],1,mean)
    plt.ylim <- range(imspe.plt)
    if("random" %in% test.types.plot){
      random.plt <- apply(performance.list$random$score[pts.important,run.vec, drop = FALSE],1,mean)
      plt.ylim <- range(c(plt.ylim, random.plt))
    }
    if("lhs" %in% test.types.plot){
      lhs.plt <- apply(performance.list$lhs$score[pts.important,run.vec, drop = FALSE],1,mean)
      plt.ylim <- range(c(plt.ylim, lhs.plt))
    }
    if("m-imspe" %in% test.types.plot){
      m.imspe.plt <- apply(performance.list$`m-imspe`$score[pts.important,run.vec, drop = FALSE],1,mean)
      plt.ylim <- range(c(plt.ylim, m.imspe.plt))
    }
    
    plot(imspe.plt, type = "l", 
         col = "blue", main = paste("average behavior after run no", i, sep = " ", collapse = " "), ylab = "score", ylim = plt.ylim)
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
    
    imspe.plt <- apply(performance.list$imspe$score[pts.important,run.vec, drop = FALSE],1,mean)
    plt.ylim <- range(imspe.plt)
    if("random" %in% test.types.plot){
      random.plt <- apply(performance.list$random$score[pts.important,run.vec, drop = FALSE],1,mean)
      plt.ylim <- range(c(plt.ylim, random.plt))
    }
    if("lhs" %in% test.types.plot){
      lhs.plt <- apply(performance.list$lhs$score[pts.important,run.vec, drop = FALSE],1,mean)
      plt.ylim <- range(c(plt.ylim, lhs.plt))
    }
    if("m-imspe" %in% test.types.plot){
      m.imspe.plt <- apply(performance.list$`m-imspe`$score[pts.important,run.vec, drop = FALSE],1,mean)
      plt.ylim <- range(c(plt.ylim, m.imspe.plt))
    }
    
    plot(imspe.plt, type = "l", 
         col = "blue", main = paste("behavior of run no", i, sep = " ", collapse = " "), ylab = "score", ylim = plt.ylim)
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
                     notes = notes, hyperparams.given = hyperparams.given, RK.vars = list(t = t, h = h), mc.iters = mc.iters, test.types = test.types,
                     candidate = candidate, uhat.starts = uhat.starts, sep = sep)
    
    
    performance.list$mcparams <- mcparams
    
    
    
    
    save(performance.list, file = file.name)
    
    
    
  }
  
}