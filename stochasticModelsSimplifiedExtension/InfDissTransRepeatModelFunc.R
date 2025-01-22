source(here::here(".//stochasticModelsSimplifiedExtension//InfDissTransModelFuncs.R"))
library(parallel)

virus_params <- function( bloodmealClearance = 1/3
                            ,propSuccessInf = 10^-3
                            ,growthRateM = 0.1
                            ,growthRateH = 0.1
                            ,carryCap = 10^8
                            ,escapeRateM = 0.05
                            ,escapeRateH = 0.05
                            ,growthVarH = v
                            ,growthVarHTrue = F
)
  return(as.list(environment()))


repeatModel <- function(x=virus_params()
    ,startingVirus=10^competenceDat$ConcMax[competenceDat$Moz %in% "Ae. aegypti"]
){

  simEachConc <- tryCatch ({ 
    simEachConc <-    lapply(startingVirus,function(conc){
      repSim <- lapply(1:100,function(y){ 
        
        if(x[9]==T){
          sc <-  as.numeric(x[8])/as.numeric(x[4])
          sh <- as.numeric(x[4])^2/ as.numeric(x[8])
          rate <- 1/sc 
          grh <- rgamma(1,shape=sh,rate=rate)
        }else{
         grh <- as.numeric(x[4]) 
        }
        
        out <- infDissTransModel(startingVirus=conc
                              ,bloodmealClearance = as.numeric(x[1])
                              ,propSuccessInf = as.numeric(x[2])
                              ,growthRateM = as.numeric(x[3])
                              ,growthRateH = grh
                              ,carryCap = as.numeric(x[5])    
                              ,escapeRateM = as.numeric(x[6])  
                              ,escapeRateH = as.numeric(x[7])
                        
                          )
        dat <- data.frame(out)
        dat$run <- y
        dat$inf <- 0
        if(dat$Mv[length(dat$Mv)]>0){dat$inf<-1}
        return(dat)
  })
      
      
      repSims <- do.call(rbind.data.frame,repSim)
      return(repSims)
      
    }) 
    
  } , error = function(e) {
    print(e)
    return(c(NA,NA,NA))
  }
  )
  
 # if(is.na(simEachConc[[1]][1])==F){  
#   simDat <- do.call(rbind.data.frame,simEachConc)
#    names(simDat) <- c("num","denom","conc")
#    return(simDat)
#  }else{
#    return(c(NA,NA,NA))
# }
}

#startTime <- Sys.time()
#test <- repeatModel(startingVirus=10^6)
#endTime <- Sys.time()
#endTime - startTime



dissSummaryFunc <- function(modelOutput){
  modelOutput$days <- round(modelOutput$time/24,1)  # create a column of .1 days
  # for each simulation (mosquito) 
  lastRunPerDay <- lapply(unique(modelOutput$run),function(y){
    runD <- modelOutput[modelOutput$run %in% y,]  # select single run
    # for each '10th of a day' in each model run select the row that has the last time point and return the results
    maxTimePerDay <- lapply(unique(runD$days),function(z){
      mtemp <- runD[runD$days %in% z,]
      maxT <- mtemp[mtemp$t %in% max(mtemp$t),]
      return(maxT[length(maxT[,1]),])
    })
    maxTimePerDay <- do.call(rbind.data.frame,maxTimePerDay) # return the results for each daily last time point
    return(maxTimePerDay)
  })
  lastRunPerDay <- do.call(rbind.data.frame,lastRunPerDay)
  # for each 10th of a day, see if any rows have Hc > 0 
  
  tempInf <- lapply(unique(lastRunPerDay$days),function(a){
    subDat <- lastRunPerDay[lastRunPerDay$days %in% a,]
    MvInf <- length(subDat$Mv[subDat$Mv >= 1])
    totalSize <- length(subDat$Mv)
    #propDiss <- HcInf/ length(unique(subDat$run))
    return(c(a,MvInf,totalSize))
  })
  
  tempDiss <- lapply(unique(lastRunPerDay$days),function(a){
    subDat <- lastRunPerDay[lastRunPerDay$days %in% a,]
    HvInf <- length(subDat$Hv[subDat$Hv >= 1])
    totalSize <- length(subDat$Hv)
    #propDiss <- HcInf/ length(unique(subDat$run))
    return(c(a,HvInf,totalSize))
  })
  tempTrans <- lapply(unique(lastRunPerDay$days),function(a){
    subDat <- lastRunPerDay[lastRunPerDay$days %in% a,]
    SvInf <- length(subDat$Sv[subDat$Sv >= 1])
    totalSize <- length(subDat$Sv)
    #propDiss <- HcInf/ length(unique(subDat$run))
    return(c(a,SvInf,totalSize))
  })
  tempInf2 <- do.call(rbind.data.frame,tempInf)
  tempDiss2 <- do.call(rbind.data.frame,tempDiss)
  tempTrans2 <- do.call(rbind.data.frame,tempTrans)
  names(tempInf2) <- c("time","numberRunsInfected","totalSize")
  names(tempDiss2) <- c("time","numberRunsDisseminated","totalSize")
  names(tempTrans2) <- c("time","numberRunsTrans","totalSize")
  dissTrans <- cbind.data.frame(tempDiss2,numberRunsTrans=tempTrans2$numberRunsTrans,numberRunsInf=tempInf2$numberRunsInfected)
  dissTrans$proportionInfected <- dissTrans$numberRunsInf/dissTrans$totalSize
  dissTrans$proportionDisseminated <- dissTrans$numberRunsDisseminated/dissTrans$totalSize
  dissTrans$proportionCapableofTransmission <- dissTrans$numberRunsTrans/dissTrans$totalSize
  return(dissTrans)
}