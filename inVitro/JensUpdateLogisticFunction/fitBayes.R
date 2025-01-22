library(R2OpenBUGS)
library(rjags)
library(coda)
library(MCMCvis)

source(here::here(".//inVitro//JensUpdateExponentialFunction//dataSecondSet.R"))
dat$logTiter <- log(dat$Titer)
dat <- dat[,-2]

ggplot(dat) +
  geom_point(aes(x=time,y=log(Titer))) +
  facet_wrap(~Chimera)

ggplot(dat) +
  geom_point(aes(x=time,y=Titer)) +
  facet_wrap(~Chimera)




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
  for(i in 1:12){                                           # for each row of data/ observation
    #mu[i] <- log(1*exp(r*time[i]))     # expected titer given r
    mu[i] <- 1*k*exp(r*time[i])/((k-1)+1*exp(r*time[i]))
    logTiter[i]~dnorm(mu[i], tau)                         # with variance
    logTiter_pred[i]~dnorm(mu[i], tau)
  }
}

# write model
model.file="model.txt"
write.model(mod,model.file)

# no initial values
inits<-list("r"=0.08,"k"=14)
# what parameters we want to track
params = c("tau","r", "k","logTiter_pred")
## hyperparameters
# number of iterations
ni = 10000
# burn in interval
nb = 1000
# thinning interval
nt = 1
# number of chains
nc = 3

c1 <- dat[dat$Chimera %in% "Chimera 1",]
c1$logTiter <- log(c1$Titer)
c1 <- c1[,c("time","logTiter")]
c1$pred <- log( (1*15*exp(0.16*c1$time))/((15-1)+1*exp(0.16*c1$time)))
c1$pred <- log(1*exp(0.16*c1$time) )

preds <- cbind.data.frame("time"=1:100,
                          "logTiter"=(1*14*exp(0.08*1:100))/((14-1)+1*exp(0.08*1:100)))

ggplot(c1) +
  geom_point(aes(x=time,y=logTiter)) +
  geom_line(data=preds,aes(x=time,y=logTiter))

c1 <- dat[dat$Chimera %in% "Chimera 1",]
c1 <- c1[,-c(1,3)]


jmod <- jags.model(file = model.file, data = c1, n.chains = nc, inits = inits, n.adapt = 1000)
#samples from the posterior for params, given MCMC hyperparameters
post = coda.samples(jmod, params, n.iter = ni, thin = nt)
# get summary of posterior samples for r parameters
MCMCsummary(post, params = c('r','k'))

# get samples from posteriors
samples = jags.samples(jmod,c('r','logTiter_pred'), length(c1$time))
# take the mean of each group of samples
posterior_means = apply(samples$logTiter_pred,1,mean)


runModel <- function(inits=NULL
                     ,params=c("tau","r","k","logTiter_pred")
                     ,ni=10000
                     ,nb=1000
                     ,nt=1
                     ,nc=3
                     ,c="Chimera 1"){
  
    chimera <- dat[dat$Chimera %in% c,]
    chimera$logTiter <- log(chimera$Titer)
    chimera <- chimera[,c("time","logTiter")]
    jmod <- jags.model(file = model.file, data = chimera, n.chains = nc, inits = inits, n.adapt = 1000)
    # iterate through jmod for the extent of the burn-in
    update(jmod, n.iter=nb, by=1)
  return(jmod)
}


c1 <- runModel(c="Chimera 1")


#****************Fit model for all chimeras**************************
modFits <- lapply(unique(dat$Chimera),function(x){
  return(runModel(c=x))
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


paramEsts <- lapply(modFits,function(x){
  return(paramEstFunc(mod=x))
})

paramEstsD <- do.call(rbind,paramEsts)
paramEstsD$param <- row.names(paramEstsD)
paramEstsD <- paramEstsD[order(paramEstsD$param),]
  
paramEstsD$chimera <- rep(unique(dat$Chimera),2)
#write.csv(paramEstsD,"draftParmsFitGCsSet2.csv")

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











