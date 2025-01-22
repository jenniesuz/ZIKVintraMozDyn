# infRate change
nllPoisInfRate <- function(loga,dat){ # Time-to-event models- exponential 
  a <- exp(loga)
  
  if (a > 0.1 | a < 0 
  ) {                                
    ll <- -9999999999999999
  }else{                                                      
    
    sim <- simPop(parms=params(probInf=a))
    
    simSub <- sim[sim$Time %in% dat$Time,]
    allDat <- merge(dat,simSub,by.x="Time",all.x=T)
    
    ll <- sum(dpois(allDat$PFU,lambda=allDat$Mv,log=T))
  }
  return(-ll)
}

#* production rate change
# infRate change
nllPoisProdRate <- function(loga,dat){ # Time-to-event models- exponential 
  a <- exp(loga)
  
  if ( a < 0 
  ) {                                
    ll <- -9999999999999999
  }else{                                                      
    
    sim <- simPop(parms=params(prodRate=a))
    
    simSub <- sim[sim$Time %in% dat$Time,]
    allDat <- merge(dat,simSub,by.x="Time",all.x=T)
    
    ll <- sum(dpois(allDat$PFU/100,lambda=allDat$Mv/100,log=T))
  }
  return(-ll)
}



# two parameter
nllPoisTwoPar <- function(loga,logb,dat){ # Time-to-event models- exponential 
  a <- exp(loga)
  b <- exp(logb)
  
  if (a > 0.1 | a < 0 | b > 50 | b < 0
      ) {                                
    ll <- -9999999999999999
  }else{                                                      
    
  sim <- simPop(parms=params(probInf=a
                             ,prodRate=b))
  
  simSub <- sim[sim$Time %in% dat$Time,]
  allDat <- merge(dat,simSub,by.x="Time",all.x=T)
  
  ll <- sum(dpois(allDat$PFU/100,lambda=allDat$Mv/100,log=T))
  }
  return(-ll)
}