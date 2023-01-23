
## Generic Runge Kutta Function

## Plug in function for differential equation

dRE.dt <- function(y = NULL,y0 = NULL, k = NULL){
  
  # k is a matrix with the columns pertaining to the reaction and the rows pertaining to the direction
  # k = [k^+_{Na}, k^+_{La}; k^-_{Na}, k^-_{La}]
  
  Na.t <- y[1]
  La.t <- y[2]

  Na.0 <- y0[1]
  La.0 <- y0[2]
  H.0 <- y0[3]
  H2A2.0 <- y0[4]
  
  LaHA.t <- La.0 - La.t
  NaHA.t <- Na.0 - Na.t
  H2A2.t <- H2A2.0 - 3*LaHA.t - NaHA.t
  H.t <- H2A2.0 + H.0 - H2A2.t
  
  dNa.dt <- k[2,1]*NaHA.t*H.t - k[1,1]*Na.t*H2A2.t
  dLa.dt <- k[2,2]*LaHA.t*H.t^3 - k[1,2]*La.t*H2A2.t^3

  return(c(dNa.dt, dLa.dt))
}

## Runge-kutta

RK <- function(y0 = NULL, k = NULL, fn = dRE.dt, t = NULL, h = NULL){
  
  steps <- floor(t/h)
  
  #y.save <<- matrix(NA, ncol = 2, nrow = steps + 1)
  
  #y0.save <<- y0
  
  y <- y0[1:2]
  
  #y.save[1,] <<- y
  
  for(i in 1:steps){
  k1 <- fn(y = y, y0 = y0, k = k)
  k2 <- fn(y = y + h*k1/2, y0 = y0, k = k)
  k3 <- fn(y = y + h*k2/2, y0 = y0, k = k)
  k4 <- fn(y = y + h*k3, y0 = y0, k = k)
  
  y <- y + (k1 +2*k2 + 2*k3 + k4)*h/6
  
  }
  
  return(y)
}

RE.conc <- function(X = NULL, U = NULL, ranges = NULL, t = NULL, h = NULL, parallel = FALSE, init.conc = NULL){
  
  X <- data.frame(X)
  names(X) <- c("NaOH.vol", "NaOH.mol","OA")
  
  elems <- c("Na", "La")
  #X <- X.code
  # if(is.data.frame(X)){
  #   X <- unname(as.matrix(X))
  # }
  # if(is.data.frame(U)){
  #   U <- unname(as.matrix(U))
  # }
  
  # Description of variables
  ## X:  Observable independent variables in coded units
  ##  Columns of X: [NaOH Volume added per unit of initial aqueous volume, NaOH mols added/L initial Aq volume, log10(O/A)]
  ## U:  Calibration parameters in coded units
  ##  Columns of U contain kinetic constants as such: [Na forward reaction, Na backward reaction, La forward reaction, La backward reaction]
  ##  Where the forward reaction is equivalent to:  La3+ + 3*H2A2 -> La(HA2)3 + 3H+
  ## ranges:  2x(ncol(X) + ncol(U)) matrix for converting coded units into natural units
  
  ### Conversion To Natural Units ###
  
  for(j in 1:ncol(X)){
    X[,j] <- (X[,j] * (ranges[2,j] - ranges[1,j])) + ranges[1,j]
  }
  for(j in 1:ncol(U)){
    U[,j] <- (U[,j] * (ranges[2,j + ncol(X)] - ranges[1,j + ncol(X)])) + ranges[1,j + ncol(X)]
  }
  
  ### Obtaining initial conditions and units required for differential equations
  
  ## Find mols of NaOH added
  
  X$OA <- 10^X$OA
  Aq.vol <- 1/(X$OA + 1)
  Org.vol <- 1-Aq.vol
  NaOH.vol <- X$NaOH.vol * Aq.vol
  
  ### Assume density of 50% NaOH is 1.52 g/cm^3
  ### Assume molarity of 50% NaOH is 18.9375 mol/L
  
  # Converts mol NaOH/Aq volume to mol NaOH
  NaOH.mol <- X$NaOH.mol*Aq.vol

  ## Adjust Na concentration to include Na and volume from NaOH addition
  
  Na.new <- (init.conc$Na*Aq.vol + NaOH.mol)/(Aq.vol + NaOH.vol)

  ## Adjust La concentration to include volume from NaOH addition
  
  La.new <- init.conc$La*Aq.vol/(Aq.vol + NaOH.vol)
  
  ## Calculate initial H+ concentration for a pH of 1.99, then neutralized by the NaOH addition.  Include NaOH volume
  
  ### working in molar basis
  H.0 <- 10^(-1.99) * Aq.vol
  pos.H <- which((H.0-NaOH.mol) > 0)
  ### working in molarity basis, including volume of NaOH addition
  H.0[pos.H] <- (H.0[pos.H] - NaOH.mol[pos.H])/(Aq.vol[pos.H] + NaOH.vol[pos.H])
  H.0[-pos.H] <- (10^(-14))/((NaOH.mol[-pos.H] - H.0[-pos.H])/(Aq.vol[-pos.H] + NaOH.vol[-pos.H]))
  
  ## Calculate H2A2 concentration in organic
  
  ## Molar Mass of EHEHPA in g/mol
  # splt <- lapply(strsplit("C16-H35-O3-P", split = "-"), strsplit, split = "")[[1]]
  # atom <- unlist(lapply(splt, function(X){X[1]}))
  # atom.number <- unlist(lapply(splt, function(X){
  #   temp <- as.numeric(paste(X[-1], sep = "", collapse = ""))
  #   if(is.na(temp)){
  #     return(1)
  #   }else{
  #     return(temp)
  #   }
  #   }))
  # ehehpa.mm <- drop(mass(atom) %*% atom.number)
  ehehpa.mm <- 306.421061
  
  ## The fresh organic is 35% EHEHPA and 65% Elixore 205 by volume
  ## Volume of ehehpa per L of fresh organic
  ## Assume:
  ### * Fresh organic is 35% ehehpa solution
  ### * ehehpa solution is 95% ehehpa
  ### * density of ehehpa solution is 0.97 g/mL
  
  molarity.ehehpa <- 0.35 * 0.95 * 0.97*1000/ehehpa.mm
  
  ## Adjust all concentrations to be on an Aq + Org volume basis, essentially adjust OA to include NaOH addition
  
  Na.0 <- Na.new*(Aq.vol + NaOH.vol)/(Aq.vol + NaOH.vol + Org.vol)
  La.0 <- La.new*(Aq.vol + NaOH.vol)/(Aq.vol + NaOH.vol + Org.vol)
  H.0 <- H.0*(Aq.vol + NaOH.vol)/(Aq.vol + NaOH.vol + Org.vol)
  H2A2.0 <- molarity.ehehpa*(Org.vol)/(Aq.vol + NaOH.vol + Org.vol)
  
  X.new <- data.frame(Na.0, La.0, H.0, H2A2.0)

  ## run RK function in parallel?

  if(nrow(X.new) == 1){

    y0 <- unname(unlist(X.new))
    k <- unname(unlist(U))
    k <- 10^k
    k <- matrix(k, ncol = 2, nrow = 2)
    y.out <- RK(y0 = y0, k = k, t = t, h = h)

    H.t <- drop(as.matrix(X.new) %*% c(1,3,1,0) - y.out %*% c(1,3))
    y.out <- c(y.out, H.t)

  }else{
    if(parallel == TRUE){
    require(doParallel)
    cores <- detectCores()
    XX <- cbind(unname(as.matrix(X.new)),unname(as.matrix(U)))
    XX <- unname(split(XX, rep(1:nrow(XX), times = ncol(XX))))
    y.out <- mclapply(X = XX, FUN = function(X, p, t, h){
      k <- X[(p +1):length(X)]
      k <- 10^k
      k <- matrix(k, ncol = 2, nrow = 2)
      y0 <- X[1:p]
      y <- RK(y0 = y0, k = k, t = t, h = h)
      return(y)
    }, p = ncol(X.new), t = t, h = h, mc.cores = cores)
    y.out <- do.call("rbind", y.out)

    H.t <- unname(as.matrix(X.new)) %*% c(1,3,1,0) - y.out %*% c(1,3)

    y.out <- cbind(y.out, H.t)
    }else{
     
      XX <- cbind(unname(as.matrix(X.new)),unname(as.matrix(U)))
      XX <- unname(split(XX, rep(1:nrow(XX), times = ncol(XX))))
      y.out <- lapply(X = XX, FUN = function(X, p, t, h){
        k <- X[(p +1):length(X)]
        k <- 10^k
        k <- matrix(k, ncol = 2, nrow = 2)
        y0 <- X[1:p]
        y <- RK(y0 = y0, k = k, t = t, h = h)
        return(y)
      }, p = ncol(X.new), t = t, h = h)
      y.out <- do.call("rbind", y.out)
      
      H.t <- unname(as.matrix(X.new)) %*% c(1,3,1,0) - y.out %*% c(1,3)
      
      y.out <- cbind(y.out, H.t)
    }
  }
  
  #provide concentrations on the basis of the total aqueous volume 
  
  y.out <- y.out * (Aq.vol + NaOH.vol + Org.vol)/(Aq.vol + NaOH.vol)
  
  
return(y.out)
  
  
}

