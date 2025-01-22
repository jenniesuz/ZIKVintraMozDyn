library(here)
source(here(".//stochasticModelsSimplified//InfDissTransRepeatModelFunc.R"))
library(grid)
library(gridExtra)
library(ggplot2)
library(parallel)
library(lhs)

# tf = 360, epilson = 0.005, starting virus = 10^6

# get combinations of probability infected and clearance that resulted in 100% midgut infection
# the original max for the probability infection is arbitrary in so much as it gives the region at which
# all simulations result in infection
midRes <- readRDS(here(".//stochasticModelsSimplified//sensitivityMidgut240917.rds"))
midRes <- midRes[midRes$proportionInfected == 1,]

# get max and min for bmClear and probabilityInfection

bmClear <- c(min(midRes$bmClear),max(midRes$bmClear))
probabilityInfection <- c(min(midRes$probabilityInfection),max(midRes$probabilityInfection))
growthRate <- c(0.0001,0.1)     # 
carryCap <- c(10^3,10^20)
escapeRate <- c(0.0000001, 0.00001)
min <- c(bmClear[1],probabilityInfection[1],growthRate[1],carryCap[1],escapeRate[1])
max <- c(bmClear[2],probabilityInfection[2],growthRate[2],carryCap[2],escapeRate[2])
params <- c("bmClear"
            ,"probabilityInfection"
            ,"growthRate"
            ,"carryCap"
            ,"escapeRate")
params <- cbind.data.frame(params,min,max)
# # select random sets of parameter values within parameter value ranges
r <- randomLHS(1000,length(params[,1]))

parmVals <- lapply(1:length(params[,1]),function(x){
  temp <- params[x,]
  randomSample <- runif(r[,x],min=temp$min,max=temp$max)
})
parmVals <- do.call(cbind.data.frame,parmVals)
names(parmVals) <- params$params

parmVals$run <- 1:1000


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
simsHS <- parLapply(cl,1:length(parmVals[,1]),function(y){    # this is lapply but in parallel - so it runs one simulation on each cluster
  
  parms <- parmVals[y,]                 # take the yth row
  
  modelOutput <- repeatModel(virus_params(bloodmealClearance = as.numeric(parms[1])
                                   ,propSuccessInf = as.numeric(parms[2])
                                   ,growthRate = as.numeric(parms[3])
                                   ,carryCap = as.numeric(parms[4])
                                   ,escapeRate = as.numeric(parms[5]))
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


saveRDS(simsHS2,here(".//stochasticModelsSimplified//sensitivityHS240917.rds"))
