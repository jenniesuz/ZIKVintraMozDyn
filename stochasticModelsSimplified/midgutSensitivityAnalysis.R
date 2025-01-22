library(here)
source(here(".//stochasticModelsSimplified//InfDissTransRepeatModelFunc.R"))
library(grid)
library(gridExtra)
library(ggplot2)
library(parallel)
library(lhs)

# tf = 168, epilson = 0.005, starting virus = 10^6

bmClear <- c(1/72,1/12)
probabilityInfection <- c(0,0.0001)
growthRate <- c(0.001,0.1)
carryCap <- c(10^3,10^20)
min <- c(bmClear[1],probabilityInfection[1],growthRate[1],carryCap[1])
max <- c(bmClear[2],probabilityInfection[2],growthRate[2],carryCap[2])
params <- c("bmClear"
            ,"probabilityInfection"
            ,"growthRate"
            ,"carryCap")
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


cl <- makeCluster(detectCores()-1)               # get cores from your computer - 1
clusterEvalQ(cl, {library(adaptivetau)})
environment(repeatModel) <- .GlobalEnv        # need to figure out whether this is needed
clusterExport(cl, varlist=c("infDissTransModel"
                            ,"virus_params"
                            ,"repeatModel"
                            ,"parmVals"
                          ),
              envir=environment())               # each cluster doesn't have access to any of the code we have run above
# so we need to send it a list of objects 


start <- Sys.time()
simsMidgut <- parLapply(cl,1:length(parmVals[,1]),function(y){    # this is lapply but in parallel - so it runs one simulation on each cluster
  
  parms <- parmVals[y,]                 # take the yth row
  
  modelOutput <- repeatModel(virus_params(bloodmealClearance = as.numeric(parms[1])
                                   ,propSuccessInf = as.numeric(parms[2])
                                   ,growthRate = as.numeric(parms[3])
                                   ,carryCap = as.numeric(parms[4])
                                   ,escapeRate = 0.00005)
                      ,startingVirus=10^6)
  
  modelOutput <- do.call(rbind,modelOutput)
  # midgut proprtion infection
  propInf <- sum(modelOutput$inf[!duplicated(modelOutput$run)])/length(modelOutput$inf[!duplicated(modelOutput$run)])
  
  return(propInf)                             # return
})
stopCluster(cl)   # stop the clusters
end <- Sys.time()
end - start

simsMidgut2 <- do.call(rbind.data.frame,simsMidgut)

simsMidgut2 <- cbind.data.frame(simsMidgut2,parmVals)
names(simsMidgut2)[1] <- "proportionInfected"

saveRDS(simsMidgut2,here(".//stochasticModelsSimplified//sensitivityMidgut240917.rds"))
