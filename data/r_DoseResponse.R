library(here)
library(ggplot2)
library(R2OpenBUGS)
library(rjags)
library(coda)
library(MCMCvis)
library(ggridges)
library(ggmcmc)

dat <- read.csv(here(".//data//30052023_Chimera_1-8_dose_response_curve_JSL.csv"))


#************Dose response curve*************
library(plyr)
library(binom)

infPrev <- ddply(dat,.(Chimera,Dose),summarise,denom=length(ID),num=sum(RT.PCR))
infPrevb <- binom.confint(x=infPrev$num,n=infPrev$denom,methods="exact")
infPrev$mean <- infPrevb$mean
infPrev$lower <- infPrevb$lower
infPrev$upper <- infPrevb$upper

infPrev$Chimera <- factor(infPrev$Chimera
                          ,labels=c("1: Senegal"
                                    ,"2: Thai"
                                    ,"3: Senegal/Thai-SP"
                                    ,"4: Thai/Senegal-SP"
                                    ,"5: Senegal/Thai-nSP"
                                    ,"6: Thai/Senegal-nSP"
                                    ,"7: Senegal/Thai-UTR"
                                    ,"8: Thai/Senegal-UTR"))



infPrev$ChimeraNSPTN <- infPrev$Chimera  
infPrev$ChimeraNSPTN <- factor(infPrev$ChimeraNSPTN
                              ,labels=c("Senegal-nSP"
                                        ,"Thai-nSP"
                                        ,"Senegal-nSP"
                                        ,"Thai-nSP"
                                        ,"Thai-nSP"
                                        ,"Senegal-nSP"
                                        ,"Senegal-nSP"
                                        ,"Thai-nSP"))



infPrev$ChimeraSPTN <- infPrev$Chimera  
infPrev$ChimeraSPTN <- factor(infPrev$ChimeraSPTN
                             ,labels=c("Senegal-SP"
                                       ,"Thai-SP"
                                       ,"Thai-SP"
                                       ,"Senegal-SP"
                                       ,"Senegal-SP"
                                       ,"Thai-SP"
                                       ,"Senegal-SP"
                                       ,"Thai-SP"))




ggplot(infPrev) +
  geom_point(aes(x=log10(Dose),y=mean,col=Chimera)) +
  geom_errorbar(aes(x=log10(Dose),ymin=lower,ymax=upper,col=Chimera)) +
  xlim(0,10) +
  ylab("Proportion of mosquitoes infected")

ggplot(infPrev) +
  geom_point(aes(x=log10(Dose),y=mean,col=ChimeraNSPTN)) +
  geom_errorbar(aes(x=log10(Dose),ymin=lower,ymax=upper,col=ChimeraNSPTN)) +
  xlim(2.5,7.5) +
  ylab("Proportion of mosquitoes infected") 



ggplot(infPrev) +
  geom_point(aes(x=log10(Dose),y=mean,col=ChimeraSPTN)) +
  geom_errorbar(aes(x=log10(Dose),ymin=lower,ymax=upper,col=ChimeraSPTN)) +
  xlim(2.5,7.5) +
  ylab("Proportion of mosquitoes infected") 



# differences in initiating infection:
# model 1: nSP is responsible
# model 2: sp is responsible (this is biologically most plausible)
# model 3: some interactive effects - e.g. Thai nsp & sp, Thai UTR & sp etc

dat$log10Dose <- log10(dat$Dose)
dat$ChimeraSPTN <- dat$Chimera  
dat$ChimeraSPTN <- factor(dat$ChimeraSPTN
                              ,labels=c("Senegal-SP"
                                        ,"Thai-SP"
                                        ,"Thai-SP"
                                        ,"Senegal-SP"
                                        ,"Senegal-SP"
                                        ,"Thai-SP"
                                        ,"Senegal-SP"
                                        ,"Thai-SP"))

test <- glm(RT.PCR~Dose+ChimeraSPTN,data=dat,family="binomial")



model= function(){
  # Likelihood
  for(i in 1:1273) {
    logit(p[i]) <- b0 + bTiter * log10Dose[i] + 
      bSP * ChimeraSPTN[i]
    RT.PCR[i] ~ dbern(p[i])
  }
  # Priors
  p0 ~ dbeta(1, 1)
  b0 <- logit(p0)
  bTiter ~ dunif(-5, 5)
  bSP ~ dunif(-5, 5)
}

# write model
model.file="doseResponseModel.txt"
write.model(model,model.file)






# no initial values
inits<-NULL
# what parameters we want to track
params = c("bSP","bTiter", "RT.PCR")
## hyperparameters
# number of iterations
ni = 10000
# burn in interval
nb = 1000
# thinning interval
nt = 1
# number of chains
nc = 3


jmod <- jags.model(file = model.file, data = dat[,c("RT.PCR","log10Dose","ChimeraSPTN")], n.chains = nc, inits = inits, n.adapt = 1000)
#samples from the posterior for params, given MCMC hyperparameters
post = coda.samples(jmod, params, n.iter = ni, thin = nt)
# get summary of posterior samples for r and k parameters
summary <- MCMCsummary(post, params = c('bSP','bTiter',"RT.PCR"))

# get samples from posteriors
samples = jags.samples(jmod,c('bSP','bTiter',"RT.PCR"), length(dat[,1]))
# take the mean of each group of samples
posterior_means = apply(samples$RT.PCR,1,mean)

rt.pcr <- summary[-c(1:2),]
rt.pcr$log10Dose <- dat$log10Dose
rt.pcr$ChimSP <- dat$ChimeraSPTN


rt.pcrS <- ddply(rt.pcr,.(ChimSP,log10Dose),summarise
                 ,denom=length(log10Dose),numM=sum(mean)
                 ,numL=sum(`2.5%`)
                 ,numU=sum(`97.5%`)
                 ,numMd=sum(`50%`))


ggplot(rt.pcrS) +
  geom_point(aes(x=log10Dose,y=numM/denom,col=ChimSP)) +
  xlim(0,10) 
