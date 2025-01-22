library(R2OpenBUGS)
library(rjags)
library(coda)
library(MCMCvis)
library(here)
library(ggridges)
library(ggmcmc)

source(here(".//inVitro//JensUpdateLogisticFunction//dataFirstSet.R"))


# Setup Model
mod = function(){
  #priors
  r~dunif(0,0.5)
  k~dunif(15,20)
  s~dunif(1,5)
  # residual variance
  tau.e<-1/(sigma.e*sigma.e)
  sigma.e~dunif(0,100) 
  # between replicate variance
  tau.u<-1/(sigma.u*sigma.u)
  sigma.u~dunif(0,100)

  #likelihood for each observational unit
  for(i in 1:12){                                                              # for each row of data/ observation
    mu[i] <- s*k*exp(r*time[i])/((k-s)+s*exp(r*time[i])) + u[Replicate[i]]     # expected titer given r
    logTiter[i]~dnorm(mu[i], tau.e)                                            # with variance
    logTiter_pred[i]~dnorm(mu[i], tau.e)
  }
  #random effects for each replicate
  for(j in 1:3){
      u[j] ~ dnorm(0,tau.u)
  }
  
}

# write model
model.file="model.txt"
write.model(mod,model.file)

# no initial values
inits<-NULL
# what parameters we want to track
params = c("tau.u","tau.e","r","k","s", "logTiter_pred")
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
c1 <- c1[,c("time","Replicate","logTiter")]

jmod <- jags.model(file = model.file, data = c1, n.chains = nc, inits = inits, n.adapt = 1000)
#samples from the posterior for params, given MCMC hyperparameters
post = coda.samples(jmod, params, n.iter = ni, thin = nt)
# get summary of posterior samples for r and k parameters
MCMCsummary(post, params = c('r','k',"s"))

# Store the chains in a data frame
postChains <- data.frame(post[[1]], iter = 1:10000)

# Use ggplot() to construct a trace plot 
ggplot(postChains, aes(x = iter, y = r)) + 
  geom_line() + theme_ridges()

ggplot(postChains, aes(x = iter, y = k)) + 
  geom_line() + theme_ridges()



# get samples from posteriors
samples = jags.samples(jmod,c('r','k','s','logTiter_pred'), length(c1$time))
# take the mean of each group of samples
posterior_means = apply(samples$logTiter_pred,1,mean)


runModel <- function(inits=list("r"=0.08,"k"=15,"s"=2.5)
                     ,params=c("tau.u","tau.e","r","k","s","logTiter_pred")
                     ,ni=10000
                     ,nb=1000
                     ,nt=1
                     ,nc=3
                     ,c="Chimera 1"){
  
    chimera <- dat[dat$Chimera %in% c,]
    chimera$logTiter <- log(chimera$Titer)
    chimera <- chimera[,c("time","Replicate","logTiter")]
    jmod <- jags.model(file = model.file, data = chimera, n.chains = nc, inits = inits, n.adapt = 1000)
    # iterate through jmod for the extent of the burn-in
    update(jmod, n.iter=nb, by=1)
  return(jmod)
}


 c1 <- runModel(c="Chimera 1")
 post = coda.samples(c1, params, n.iter = ni, thin = nt)
 # get summary of posterior samples for r parameters
 summary <- MCMCsummary(post, params = c('r','k','s'))
# 
# 1*k*exp(r*time[i])/((k-1)+1*exp(r*time[i]))
# 
 mean <- sapply(1:100
                ,FUN=function(x){summary$mean[3]*summary$mean[2]*exp(summary$mean[1]*x)/((summary$mean[2]-summary$mean[3])+summary$mean[3]*exp(summary$mean[1]*x))}
 )
 lower <- sapply(1:100
                 ,FUN=function(x){summary$mean[3]*summary$`2.5%`[2]*exp(summary$`2.5%`[1]*x)/((summary$`2.5%`[2]-summary$mean[3])+summary$mean[3]*exp(summary$`2.5%`[1]*x))}
 )
 upper <- sapply(1:100
                 ,FUN=function(x){summary$mean[3]*summary$`97.5%`[2]*exp(summary$`97.5%`[1]*x)/((summary$`97.5%`[2]-summary$mean[3])+summary$mean[3]*exp(summary$`97.5%`[1]*x))}
 )
# 
 predDat <- cbind.data.frame(time=1:100,mean,lower,upper)
 ggplot() +
   geom_point(data=dat[dat$Chimera %in% "Chimera 1",],aes(x=time,y=log(Titer))) +
   geom_line(data=predDat,aes(x=time,y=mean),col="blue")  +
   geom_ribbon(data=predDat,aes(x=time,ymin =lower, ymax =upper),fill="blue", alpha = .2) 
# 


#****************Fit model for all chimeras**************************
modFits <- lapply(unique(dat$Chimera),function(x){
  return(runModel(c=x))
})
names(modFits) <- unique(dat$Chimera)
#****************************************************************



postSamples <- lapply(1:length(modFits),function(x){
  temp <- modFits[[x]]
  samples <- coda.samples(temp, variable.names=c("r","k","s"), n.iter=1000, thin = 1)
  samples <- do.call(rbind.data.frame,samples)
  samples$Chimera <- as.character(names(modFits)[x])
  return(samples)
})

postChains <- lapply(1:length(modFits),function(x){
  temp <- modFits[[x]]
  samples <- coda.samples(temp, variable.names=c("r","k","s"), n.iter=10000, thin = 1)
  samples <- data.frame(samples[[1]], iter = 1:10000)
  samples$Chimera <- as.character(names(modFits)[x])
  return(samples)
})

postSamples <- do.call(rbind.data.frame,postSamples)

postChains <- do.call(rbind.data.frame,postChains)
write.csv(postChains,here(".\\inVitro\\JensUpdateLogisticFunction\\postDist1stSetConstantGrowthRateCh.csv"))

write.csv(postSamples,here(".\\inVitro\\JensUpdateLogisticFunction\\postDist1stSetConstantGrowthRate.csv"))

#*************************Plot posterior param*********************
postSamples$Chimera <- factor(postSamples$Chimera
                              ,labels=c("1: Senegal"
                                        ,"2: Thai"
                                        ,"3: Senegal/Thai-SP"
                                        ,"4: Thai/Senegal-SP"
                                        ,"5: Senegal/Thai-nSP"
                                        ,"6: Thai/Senegal-nSP"
                                        ,"7: Senegal/Thai-UTR"
                                        ,"8: Thai/Senegal-UTR"))




#tiff(here(".//inVitro//JensUpdateLogisticFunction//fig_1stSetEstimateR.tiff"), height =5, width = 4, units = 'in', compression="lzw", res=400)
ggplot(data = postSamples, aes(x = r, y = Chimera)) + 
  stat_density_ridges(quantile_lines = TRUE, 
                      quantiles = c(0.025, 0.5, 0.975), 
                      alpha = 0.2
                      ,fill="blue") + 
  labs(title="Growth rate (r)") +
  xlab("Posterior estimate") + 
  ylab("") +
  theme_bw() +
  theme(panel.border = element_blank()                   
        ,axis.line = element_line(color = 'black')
        ,text=element_text(size=10)
        ,plot.margin=unit(c(0.2,0.1,0.1,0.1), "cm")
        ,axis.text=element_text(size=8)
        ,legend.key.size = unit(0.5,"line")  
        ,legend.background = element_blank()
        ,legend.text=element_text(size=8)
        ,strip.text=element_text(size=8)
        ,strip.background =element_blank()
  )
#dev.off()

#tiff(here(".//inVitro//JensUpdateLogisticFunction//fig_1stSetEstimatek.tiff"), height =5, width = 4, units = 'in', compression="lzw", res=400)
ggplot(data = postSamples, aes(x = k, y = Chimera)) + 
  stat_density_ridges(quantile_lines = TRUE, 
                      quantiles = c(0.025, 0.5, 0.975), 
                      alpha = 0.2
                      ,fill="blue") + 
  labs(title="Carrying capacity (k)") +
  xlab("Posterior estimate") + 
  ylab("") +
  theme_bw() +
  theme(panel.border = element_blank()                   
        ,axis.line = element_line(color = 'black')
        ,text=element_text(size=10)
        ,plot.margin=unit(c(0.2,0.1,0.1,0.1), "cm")
        ,axis.text=element_text(size=8)
        ,legend.key.size = unit(0.5,"line")  
        ,legend.background = element_blank()
        ,legend.text=element_text(size=8)
        ,strip.text=element_text(size=8)
        ,strip.background =element_blank()
  )
#dev.off()


#tiff(here(".//inVitro//JensUpdateLogisticFunction//fig_1stSetEstimateS.tiff"), height =5, width = 4, units = 'in', compression="lzw", res=400)
ggplot(data = postSamples, aes(x = s, y = Chimera)) + 
  stat_density_ridges(quantile_lines = TRUE, 
                      quantiles = c(0.025, 0.5, 0.975), 
                      alpha = 0.2
                      ,fill="blue") + 
  labs(title="Starting virus (s)") +
  xlab("Posterior estimate") + 
  ylab("") +
  theme_bw() +
  theme(panel.border = element_blank()                   
        ,axis.line = element_line(color = 'black')
        ,text=element_text(size=10)
        ,plot.margin=unit(c(0.2,0.1,0.1,0.1), "cm")
        ,axis.text=element_text(size=8)
        ,legend.key.size = unit(0.5,"line")  
        ,legend.background = element_blank()
        ,legend.text=element_text(size=8)
        ,strip.text=element_text(size=8)
        ,strip.background =element_blank()
  )
#dev.off()





#*************get rate estimate and credible intervals*************
paramEstFunc <- function(mod){
  #samples from the posterior for params, given MCMC hyperparameters
  post = coda.samples(mod, params, n.iter = ni, thin = nt)
  # get summary of posterior samples for r parameters
  summary <- MCMCsummary(post, params = c('r','k','s'))
  return(summary)
}


paramEsts <- lapply(modFits,function(x){
  return(paramEstFunc(mod=x))
})

paramEstsD <- do.call(rbind,paramEsts)

paramEstsD$param <- rep(c("r","k","s"),length(unique(dat$Chimera)))

paramEstsD <- paramEstsD[order(paramEstsD$param),]

paramEstsD$chimera <- rep(unique(dat$Chimera),3)



#***************Predictions************************
predsFunc <- function(Chimera){
  summary <- paramEstsD[paramEstsD$chimera %in% Chimera,]
  mean <- sapply(1:100
                 ,FUN=function(x){summary$mean[3]*summary$mean[1]*exp(summary$mean[2]*x)/((summary$mean[1]-summary$mean[3])+summary$mean[3]*exp(summary$mean[2]*x))}
  )
  lower <- sapply(1:100
                  ,FUN=function(x){summary$`2.5%`[3]*summary$`2.5%`[1]*exp(summary$`2.5%`[2]*x)/((summary$`2.5%`[1]-summary$`2.5%`[3])+summary$`2.5%`[3]*exp(summary$`2.5%`[2]*x))}
  )
  upper <- sapply(1:100
                  ,FUN=function(x){summary$`97.5%`[3]*summary$`97.5%`[1]*exp(summary$`97.5%`[2]*x)/((summary$`97.5%`[1]-summary$`97.5%`[3])+summary$`97.5%`[3]*exp(summary$`97.5%`[2]*x))}
  )
  
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

predData$Chimera <-  factor(predData$Chimera
                            ,labels=c("1: Senegal"
                                      ,"2: Thai"
                                      ,"3: Senegal/Thai-SP"
                                      ,"4: Thai/Senegal-SP"
                                      ,"5: Senegal/Thai-nSP"
                                      ,"6: Thai/Senegal-nSP"
                                      ,"7: Senegal/Thai-UTR"
                                      ,"8: Thai/Senegal-UTR"))

dat$Chimera <-  factor(dat$Chimera
                       ,labels=c("1: Senegal"
                                 ,"2: Thai"
                                 ,"3: Senegal/Thai-SP"
                                 ,"4: Thai/Senegal-SP"
                                 ,"5: Senegal/Thai-nSP"
                                 ,"6: Thai/Senegal-nSP"
                                 ,"7: Senegal/Thai-UTR"
                                 ,"8: Thai/Senegal-UTR"))



dat$ChimeraPTN <- dat$Chimera  
dat$ChimeraPTN <- factor(dat$ChimeraPTN
                         ,labels=c("Senegal-nSP"
                                   ,"Thai-nSP"
                                   ,"Senegal-nSP"
                                   ,"Thai-nSP"
                                   ,"Thai-nSP"
                                   ,"Senegal-nSP"
                                   ,"Senegal-nSP"
                                   ,"Thai-nSP"))

predData$ChimeraPTN <- predData$Chimera  
predData$ChimeraPTN <- factor(predData$ChimeraPTN
                         ,labels=c("Senegal-nSP"
                                   ,"Thai-nSP"
                                   ,"Senegal-nSP"
                                   ,"Thai-nSP"
                                   ,"Thai-nSP"
                                   ,"Senegal-nSP"
                                   ,"Senegal-nSP"
                                   ,"Thai-nSP"))



#tiff(here(".//inVitro//JensUpdateLogisticFunction//fig_1stSetModelFits.tiff"), height =5, width = 5, units = 'in', compression="lzw", res=400)
ggplot() +
  geom_point(data=dat,aes(x=time,y=log(Titer))) +
  geom_line(data=predData,aes(x=time,y=mean),col="blue")  +
  geom_ribbon(data=predData,aes(x=time,ymin =lower, ymax =upper),fill="blue", alpha = .2) +
  facet_wrap(~Chimera) +
  theme_bw() +
  theme(panel.border = element_blank()                   
        ,axis.line = element_line(color = 'black')
        ,text=element_text(size=10)
        ,plot.margin=unit(c(0.2,0.1,0.1,0.1), "cm")
        ,axis.text=element_text(size=8)
        ,legend.key.size = unit(0.5,"line")  
        ,legend.background = element_blank()
        ,legend.text=element_text(size=8)
        ,strip.text=element_text(size=8)
        ,strip.background =element_blank()
  )
#dev.off()


#tiff(here(".//inVitro//JensUpdateLogisticFunction//fig_1stSetModelFitsshowingbyNSP.tiff"), height =4, width = 5, units = 'in', compression="lzw", res=400)
ggplot() +
  geom_point(data=dat,aes(x=time,y=log(Titer),col=Chimera)) +
  geom_line(data=predData,aes(x=time,y=mean,col=Chimera),size=1)  +
  geom_ribbon(data=predData,aes(x=time,ymin =lower, ymax =upper, fill=Chimera), alpha = .05) +
  facet_wrap(~ChimeraPTN) +
  theme_bw() +
  theme(panel.border = element_blank()                   
        ,axis.line = element_line(color = 'black')
        ,text=element_text(size=10)
        ,plot.margin=unit(c(0.2,0.1,0.1,0.1), "cm")
        ,axis.text=element_text(size=8)
        ,legend.key.size = unit(0.5,"line")  
        ,legend.background = element_blank()
        ,legend.text=element_text(size=8)
        ,strip.text=element_text(size=8)
        ,strip.background =element_blank()
  )
#dev.off()


