
nll <- function(loga,logb,logc,logd,dat){ # Time-to-event models- exponential 

    sim <- simPop(parms=params(infRate=exp(loga)
                               ,muV=exp(logb)
                               ,scalingParameter=exp(logc)))
    
    simSub <- sim[sim$time %in% dat$time,]
    allDat <- merge(dat,simSub,by.x="time",all.x=T)
    
    ll <- sum(dnorm(log(allDat$Titer),mean=log(allDat$V),sd=exp(logd),log=T))

  return(-ll)
}




nll2 <- function(loga,logc,logd,dat){ # Time-to-event models- exponential 
  
  sim <- simPop(parms=params(infRate=exp(loga)
                             ,scalingParameter=exp(logc)))
  
  simSub <- sim[sim$time %in% dat$time,]
  allDat <- merge(dat,simSub,by.x="time",all.x=T)
  
  ll <- sum(dnorm(log(allDat$Titer),mean=log(allDat$V),sd=exp(logd),log=T))
  
  return(-ll)
}

