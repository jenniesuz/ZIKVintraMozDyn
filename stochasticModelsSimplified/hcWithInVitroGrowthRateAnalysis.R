
# in the laboratory results all chimeras produced a disseminated infection across
# all time points but transmission looked different 

library(here)
source(here(".//stochasticModelsSimplified//InfDissTransRepeatModelFunc.R"))
library(grid)
library(gridExtra)
library(ggplot2)
library(parallel)
library(lhs)



#******************Set up clusters*********************
cl <- makeCluster(detectCores()-1)               # get cores from your computer - 1
clusterEvalQ(cl, {library(adaptivetau)})
environment(repeatModel) <- .GlobalEnv        # need to figure out whether this is needed
clusterExport(cl, varlist=c("infDissTransModel"
                            ,"virus_params"
                            ,"repeatModel"
                            ,"parmVals"
                            ,"dissSummaryFunc"
),
envir=environment())               # each cluster doesn't have access to any of the code we have run above
# so we need to send it a list of objects 

#*****************Run simulations**********************
start <- Sys.time()
simsHS <- parLapply(cl,1:length(parmVals),function(y){    # this is lapply but in parallel - so it runs one simulation on each cluster
  
  parms <- parmVals[y]                 # take the yth row
  
  modelOutput <- repeatModel(virus_params(bloodmealClearance = 1/72
                                          ,propSuccessInf = 0.001
                                          ,growthRate = 0.04
                                          ,carryCap = 10^20
                                          ,escapeRate = 0.0005) # 0.001
                             ,startingVirus=10^6)
  
  modelOutput <- do.call(rbind,modelOutput)
  # get summaries of proportion disseminated and capable of transmitting at the 
  # same time points the experiment was carried out - 7, 10 14 days
  props <- dissSummaryFunc(modelOutput)
  props$roundTime <- round(props$time,0)
  sevenDaysH <- max(props$proportionDisseminated[props$roundTime == 7])
  tenDaysH <- max(props$proportionDisseminated[props$roundTime == 10])
  fourteenDaysH <- max(props$proportionDisseminated[props$roundTime == 14])
  sevenDaysS <- max(props$proportionCapableofTransmission[props$roundTime == 7])
  tenDaysS <- max(props$proportionCapableofTransmission[props$roundTime == 10])
  fourteenDaysS <- max(props$proportionCapableofTransmission[props$roundTime == 14])
  
  propsHC <- c(sevenDaysH,tenDaysH,fourteenDaysH,sevenDaysS,tenDaysS,fourteenDaysS)
  return(propsHC)                             # return
})
stopCluster(cl)   # stop the clusters
end <- Sys.time()
end - start


#***************Sort and save output*****************

simsHS2 <- do.call(rbind.data.frame,simsHS)

simsHS2 <- cbind.data.frame(simsHS2,parmVals)


saveRDS(simsHS2,"HSSenegalInvitroGrowthRatees0.0005.rds")




#***************************************THAI***************************************************
# use sensitivity analysis to guide estimates of fixed parameters - want max escape rate for now
# as lab results dissemination had occurred already at 7 days
parmVals <- postDist$r[postDist$Chimera %in% "2: Thai"]
parmVals <- sample(parmVals,1000,replace=F)
hist(parmVals)

#******************Set up clusters*********************
cl <- makeCluster(detectCores()-1)               # get cores from your computer - 1
clusterEvalQ(cl, {library(adaptivetau)})
environment(repeatModel) <- .GlobalEnv        # need to figure out whether this is needed
clusterExport(cl, varlist=c("infDissTransModel"
                            ,"virus_params"
                            ,"repeatModel"
                            ,"parmVals"
                            ,"dissSummaryFunc"
),
envir=environment())               # each cluster doesn't have access to any of the code we have run above
# so we need to send it a list of objects 

#*****************Run simulations**********************
start <- Sys.time()
simsHS <- parLapply(cl,1:length(parmVals),function(y){    # this is lapply but in parallel - so it runs one simulation on each cluster
  
  parms <- parmVals[y]                 # take the yth row
  
  modelOutput <- repeatModel(virus_params(bloodmealClearance = 1/72
                                          ,propSuccessInf = 0.001
                                          ,growthRate = as.numeric(parms)
                                          ,carryCap = 10^20
                                          ,escapeRate = 0.0005) #0.001
                             ,startingVirus=10^6)
  
  modelOutput <- do.call(rbind,modelOutput)
  # get summaries of proportion disseminated and capable of transmitting at the 
  # same time points the experiment was carried out - 7, 10 14 days
  props <- dissSummaryFunc(modelOutput)
  props$roundTime <- round(props$time,0)
  sevenDaysH <- max(props$proportionDisseminated[props$roundTime == 7])
  tenDaysH <- max(props$proportionDisseminated[props$roundTime == 10])
  fourteenDaysH <- max(props$proportionDisseminated[props$roundTime == 14])
  sevenDaysS <- max(props$proportionCapableofTransmission[props$roundTime == 7])
  tenDaysS <- max(props$proportionCapableofTransmission[props$roundTime == 10])
  fourteenDaysS <- max(props$proportionCapableofTransmission[props$roundTime == 14])
  
  propsHC <- c(sevenDaysH,tenDaysH,fourteenDaysH,sevenDaysS,tenDaysS,fourteenDaysS)
  return(propsHC)                             # return
})
stopCluster(cl)   # stop the clusters
end <- Sys.time()
end - start


#***************Sort and save output*****************

simsHS2 <- do.call(rbind.data.frame,simsHS)

simsHS2 <- cbind.data.frame(simsHS2,parmVals)


saveRDS(simsHS2,"HSThaiInvitroGrowthRatees0.0005.rds")
