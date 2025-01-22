library(R2OpenBUGS)
library(rjags)
library(coda)
library(MCMCvis)

# Idea would be to fit the individual replicates then do Kruskall-Wallis to compare rates
# Unless Mike has better ideas - will ask for advice

# same as fitBayes script but now for individual reps
source(here::here(".//inVitro//JensUpdateLogisticFunction//dataSecondSet.R"))
dat$logTiter <- log(dat$Titer)
dat$ChimeraRep <- paste(dat$Chimera,dat$Replicate,sep=" ")

# Setup Model
mod = function(){
  #priors
  #r~dnorm(0,.001)        # growth rate
  r~dunif(0,0.5)
  k~dunif(10,15)
  # residual variance
  tau<-1/(sigma*sigma)
  sigma~dunif(0,100) 

  #likelihood for each observational unit
  for(i in 1:4){                                           # for each row of data/ observation
    #mu[i] <- log(1*exp(r*time[i]))     # expected titer given r
    mu[i] <- 1*k*exp(r*time[i])/((k-1)+1*exp(r*time[i]))
    logTiter[i]~dnorm(mu[i], tau)                         # with variance
    logTiter_pred[i]~dnorm(mu[i], tau)
  }
}

# write model
model.file="model.txt"
write.model(mod,model.file)

#************function to loop over chimeras and replicates***********
runModel <- function(inits=NULL
                     ,params=c("tau","r","k","logTiter_pred")
                     ,ni=10000
                     ,nb=1000
                     ,nt=1
                     ,nc=3
                     ,chimera){
  
    chimera$logTiter <- log(chimera$Titer)
    chimera <- chimera[,c("time","logTiter")]
    jmod <- jags.model(file = model.file, data = chimera, n.chains = nc, inits = inits, n.adapt = 1000)
    # iterate through jmod for the extent of the burn-in
    update(jmod, n.iter=nb, by=1)
  return(jmod)
}



#****************Fit model for all chimeras**************************
modFits <- lapply(unique(dat$ChimeraRep),function(x){
  data <- dat[dat$ChimeraRep %in% x,]
  return(runModel(chimera=data))
})

#****************************************************************


#*************get rate estimate and credible intervals*************
paramEstFunc <- function(mod){
  #samples from the posterior for params, given MCMC hyperparameters
  post = coda.samples(mod, params, n.iter = ni, thin = nt)
  # get summary of posterior samples for r parameters
  summary <- MCMCsummary(post, params = c('r',"k"))
  return(summary)
}

chimrep <- lapply(unique(dat$ChimeraRep),function(x){
  return(rep(x,2))
})
chimrep <- unlist(chimrep)

chim <- lapply(unique(dat$Chimera),function(x){
  return(rep(x,6))
})
chim <- unlist(chim)

paramEsts <- lapply(modFits,function(x){
  return(paramEstFunc(mod=x))
})

paramEstsD <- do.call(rbind,paramEsts)
paramEstsD$param <- row.names(paramEstsD)
paramEstsD$chimRep <- chimrep
paramEstsD$chimera <- chim  
#write.csv(paramEstsD,"draftParmsFitGCsSet2IndReps.csv")



#**********or would you compare the posterior distributions in some way instead??****
paramEstsD$param <- rep(c("r","k"),27)
test <- kruskal.test(paramEstsD$mean[paramEstsD$param %in% "r"],g=paramEstsD$chimera[paramEstsD$param %in% "r"])



#***************Predictions************************
predsFunc <- function(Chimera){
  summary <- paramEstsD[paramEstsD$chimera %in% Chimera,]
  mean <- sapply(1:100
                 ,FUN=function(x){1*summary$mean[1]*exp(summary$mean[2]*x)/((summary$mean[1]-1)+1*exp(summary$mean[2]*x))}
                   )
  lower <- sapply(1:100
                 ,FUN=function(x){1*summary$`2.5%`[1]*exp(summary$`2.5%`[2]*x)/((summary$`2.5%`[1]-1)+1*exp(summary$`2.5%`[2]*x))}
  )
  upper <- sapply(1:100
                 ,FUN=function(x){1*summary$`97.5%`[1]*exp(summary$`97.5%`[2]*x)/((summary$`97.5%`[1]-1)+1*exp(summary$`97.5%`[2]*x))}
  )
  
  #mean <- sapply(1:100, FUN=function(x){1*exp(summary$mean*x)})
  #lower <- sapply(1:100, FUN = function(x){1*exp(summary$`2.5%`*x)})
  #upper <- sapply(1:100, FUN = function(x){1*exp(summary$`97.5%`*x)})
  predDat <- cbind.data.frame(time=1:100,mean,lower,upper)
  return(predDat)
}

preds <- lapply(unique(paramEstsD$chimera),function(x){
  return(predsFunc(Chimera=x))
})


chims <- lapply(unique(dat$Chimera),function(x){
  return(rep(x,100))
})

preds <- do.call(rbind.data.frame,preds)
chims <- unlist(chims)
predData <- cbind.data.frame(preds,chims)
names(predData)[5] <- "Chimera"

ggplot() +
  geom_point(data=dat,aes(x=time,y=log(Titer))) +
  geom_line(data=predData,aes(x=time,y=mean),col="blue")  +
 geom_ribbon(data=predData,aes(x=time,ymin =lower, ymax =upper),fill="blue", alpha = .2) +
  facet_wrap(~Chimera)


ggplot() +
  geom_point(data=dat,aes(x=time,y=log(Titer),col=Chimera)) +
  geom_line(data=predData,aes(x=time,y=mean,col=Chimera,group=Chimera))  +
  geom_ribbon(data=predData,aes(x=time,ymin =lower, ymax =upper,fill=Chimera,group=Chimera), alpha = .2) 











