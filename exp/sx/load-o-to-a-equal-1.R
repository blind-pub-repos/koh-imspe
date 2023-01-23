#########
### SCRIPT FOR LOADING DATA WITH O/A OF 1
#########

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

eps <- sqrt(.Machine$double.eps)

#############################################################################
### BEGIN Field data organization
#############################################################################


dat <- read.csv("OA-1-to-1.csv")
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


### Delete test where it was noted there was a longer mixing time

dat <- dat[-((which(dat$Sample == "K001-220829KA2")-1):which(dat$Sample == "K001-220829KA2")),]

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

### Add column for vol of NaOH mols added/L of Aqueous and column for NaOH volume used/L of Aqueous

### Add column for vol and mol of NaOH used
### Assume 
####  density of 50% NaOH is 1.52 g/cm^3
####  molarity of 50% NaOH is 18.9375 mol/L
####  density of 5% NaOH is 1.0554 g/cm^3
####  molarity of 5% NaOH is 1.32 mol/L

dat$X50..NaOH.Used[is.na(dat$X50..NaOH.Used)] <- 0

dat$NaOH.vol <- NA
dat$NaOH.mol <- NA

dat$NaOH.vol[dat$Sample.Type == "Feed"] <- (dat$X5..NaOH.Used[dat$Sample.Type == "Aqueous"]/1.0554 + 
                                              dat$X50..NaOH.Used[dat$Sample.Type == "Aqueous"]/1.52)/dat$Aq.In[dat$Sample.Type == "Aqueous"]

dat$NaOH.mol[dat$Sample.Type == "Feed"] <- (dat$X5..NaOH.Used[dat$Sample.Type == "Aqueous"]/1.0554*1.32 + 
                                              dat$X50..NaOH.Used[dat$Sample.Type == "Aqueous"]/1.52*18.9375)/dat$Aq.In[dat$Sample.Type == "Aqueous"]




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

X.1to1 <- X
Y.1to1 <- Y