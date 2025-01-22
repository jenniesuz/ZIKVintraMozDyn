library(R2OpenBUGS)
library(rjags)
library(coda)
library(MCMCvis)
library(here)
library(ggridges)
library(ggmcmc)
library(grid)
library(gridExtra)
library(BayesTools)
library(runjags)
library(parallel)

source(here(".//inVitro//JensUpdateLogisticFunction//dataFirstSet.R"))

# Depaoli et al (2020) the importance of prior sensitivity analysis in Bayesian
# statistics: demonstrations using an interactive shiny app

# one of the biggest aids for examining the role of priors
# can be to visually examine the resulting posterior distributions
# across different prior settings. 
# visual aids are particularly important
# a separate sensitivity analysis can be conducted 
# on each parameter and another examines the combined
# specification

# r range 0, 0.5
# k range 15, 20
# s range 1, 5

# assume gamma priors

#*********************starting virus********************
m <- 3
v <- 5
sc <- v/m
sh <- m^2/v
rate <- 1/sc
x.max = qgamma(0.999, shape=sh, scale=sc)
x = seq(from=0, to=x.max, by=x.max/1000)
dens = dgamma(x, shape=sh, scale=sc)
plot(x, dens, type='l')
# dividing variance by mean gives scale parameter
# shape = mean*scale
# rate = 1/scale

# mean 3, var 5, sh 1.8, rate 0.6


#******************carrying capacity**************************
m <- 15
v <- 40  
sc <- v/m
sh <- m^2/v
rate <- 1/sc
x.max = qgamma(0.999, shape=sh, scale=sc)
x = seq(from=0, to=x.max, by=x.max/1000)
dens = dgamma(x, shape=sh, scale=sc)
plot(x, dens, type='l')
# dividing variance by mean gives scale parameter
# shape = mean*scale
# rate = 1/scale

# sh 11.25, rate 0.75


#********************************growth rate********************
m <- 0.25
v <- 0.1
sc <- v/m
sh <- m^2/v
rate <- 1/sc
x.max = qgamma(0.999, shape=sh, scale=sc)
x = seq(from=0, to=x.max, by=x.max/1000)
dens = dgamma(x, shape=sh, scale=sc)
plot(x, dens, type='l')

# mean 0.25, var 0.05, sh 1.25 rate 5


#*************************************************************

# halving and doubling parameters
startingVMeans <- c(3,6,1.5)
startingVvar <- c(5,10,2.5)

ccMeans <- c(15,30,7.5)
ccVar <- c(20,40,10)

grMeans <- c(0.25,0.5,0.125)
grVar <- c(0.05,0.1,0.025)

pcombs <- expand.grid(startingVMeans,startingVvar,ccMeans,ccVar,grMeans,grVar)
names(pcombs) <- c("startingVMeans","startingVvar","ccMeans","ccVar","grMeans","grVar")

pcombs <- pcombs[1:10,]

# define likelihood for the data
model_syntax <-
  "model{
      tau.e<-1/(sigma.e*sigma.e)
      sigma.e~dunif(0,100) 
      
      tau.u<-1/(sigma.u*sigma.u)
      sigma.u~dunif(0,100)

     
    for(i in 1:12){                                                             
      mu[i] <- s*k*exp(r*time[i])/((k-s)+s*exp(r*time[i])) + u[Replicate[i]]     
      logTiter[i]~dnorm(mu[i], tau.e)                                           
      logTiter_pred[i]~dnorm(mu[i], tau.e)
    }
   
    for(j in 1:3){
      u[j] ~ dnorm(0,tau.u)
    }
    }"

#************************data***********************************

dat$logTiter <- log(dat$Titer)
chimera <- dat[dat$Chimera %in% "Chimera 1",]
chimera <- dat[,c("time","Replicate","logTiter")]

#*****************Original priors******************************
#***************************Using BayesTools******************************


p0r <- prior(distribution = "gamma",  parameters = list(shape=1.25, rate=5))
p0s <- prior(distribution = "gamma", parameters=list(shape=1.8,rate=0.6))
p0k <- prior(distribution = "gamma", parameters=list(shape=11.25,rate=0.75))

psige <- prior(distribution = "uniform", parameters=list(a=0,b=2))
psigu <- prior(distribution = "uniform", parameters=list(a=0,b=2))

priors_list0 <- list(r = p0r, s = p0s, k = p0k, sigma.e = psige, sigma.u = psigu)

# define likelihood for the data
model_syntax <-
  "model{
      tau.e<-1/(sigma.e*sigma.e)
    
      tau.u<-1/(sigma.u*sigma.u)

    for(i in 1:12){                                                             
      mu[i] <- s*k*exp(r*time[i])/((k-s)+s*exp(r*time[i])) + u[Replicate[i]]     
      logTiter[i]~dnorm(mu[i], tau.e)                                           
      logTiter_pred[i]~dnorm(mu[i], tau.e)
    }
   
    for(j in 1:3){
      u[j] ~ dnorm(0,tau.u)
    }
    }"


# fit the models
fit0 <- JAGS_fit(model_syntax, data=chimera, priors_list0, seed = 0)

estFit0 <- runjags_estimates_table(fit0)

# define log posterior for bridge sampling
log_posterior <- function(parameters=data.frame("r"=mean(p0r)
                                                ,"s"=mean(p0s)
                                                ,"k"=mean(p0k)
                                                ,"sigma.e"=mean(psige)
                                                ,"sigma.u"=mean(psigu))
                          ,data=chimera){
  
  tau.e<-1/(parameters$sigma.e*parameters$sigma.e)
  tau.u<-1/(parameters$sigma.u*parameters$sigma.u)
  
  mu <- rep(NA,12)
  logTiter_pred <- rep(NA,12)
  u <- rep(NA,3)
  for(j in 1:3){
    u[j] <- dnorm(1,mean=0,sd=tau.u)
  }
  
  for(i in 1:12){                                                             
    mu[i] <- parameters$s*parameters$k*exp(parameters$r*data$time[i])/((parameters$k-parameters$s)+parameters$s*exp(parameters$r*data$time[i])) + u[data$Replicate[i]]     
    logTiter_pred[i] <- dnorm(mu[i], tau.e)
  }
  
 sum(dnorm(data$logTiter, logTiter_pred, 1, log = TRUE))
}


log_posterior()

prior_prob <- function(parameters=data.frame("r"=mean(p0r)
                                             ,"s"=mean(p0s)
                                             ,"k"=mean(p0k)
                                             ,"sigma.e"=mean(psige)
                                             ,"sigma.u"=mean(psigu))
                       ,data=chimera){
  tau.e<-1/(parameters$sigma.e*parameters$sigma.e)
  tau.u<-1/(parameters$sigma.u*parameters$sigma.u)
  
  mu <- rep(NA,12)
  logTiter_pred <- rep(NA,12)
  u <- rep(NA,3)
  for(j in 1:3){
    u[j] <- dnorm(1,mean=0,sd=tau.u)
  }
  
  for(i in 1:12){                                                             
    mu[i] <- parameters$s*parameters$k*exp(parameters$r*data$time[i])/((parameters$k-parameters$s)+parameters$s*exp(parameters$r*data$time[i])) + u[data$Replicate[i]]     
    logTiter_pred[i] <- dnorm(mu[i], tau.e)
  }
  
  sum(dnorm(data$logTiter, logTiter_pred, 1, log = F))
  
}

prior_prob()

# get marginal likelihoods
marglik0 <- list(
  logml = log_posterior()
)

class(marglik0) <- "bridge"
marglik0 <- JAGS_bridgesampling(fit=fit0, data=chimera, prior_list=priors_list0, log_posterior)
prob0 <- prior_prob()

pcombs <- pcombs[1:5,]
#***********************Sensitivity****************
cl <- makeCluster(detectCores()-1)               # get cores from your computer - 1

environment(log_posterior) <- .GlobalEnv        
clusterExport(cl, varlist=c("pcombs"
                            ,"chimera"
                            ,"prior"
                            ,"model_syntax"
                            ,"psige"
                            ,"psigu"
                            ,"marglik0"
                            ,"prob0"
                            ,"prior_prob"
                            ,"log_posterior"
                            ,"JAGS_fit"
                            ,"runjags_estimates_table"
                            ,"JAGS_bridgesampling"
                            ,"inclusion_BF"
                           ),
              envir=environment())


start <- Sys.time()
sims <- parLapply(cl,1:length(pcombs[,1]),function(y){    # this is lapply but in parallel - so it runs one simulation on each cluster
  y <- as.numeric(pcombs[y,])
  sSh <- y[1]^2/y[2]
  ssc <- y[2]/y[1]
  sRate <- 1/ssc
  
  kSh <- y[3]^2/y[4]
  ksc <- y[4]/y[3]
  kRate <- 1/ksc
  
  rSh <- y[5]^2/y[6]
  rsc <- y[6]/y[5]
  rRate <- 1/rsc
  
  #*********************************************************************
  priorr <- prior(distribution="gamma",parameters=list(shape=rSh,rate=rRate))
  priors <- prior(distribution="gamma",parameters=list(shape=sSh,rate=sRate))
  priork <- prior(distribution="gamma",parameters=list(shape=kSh,rate=kRate))
  
  priors_list <- list(r = priorr, s = priors, k = priork, sigma.e = psige, sigma.u = psigu)
  
  # fit the models
  fit1 <- JAGS_fit(model_syntax, data=chimera, priors_list, seed = 0)
  estFit1 <- runjags_estimates_table(fit1)
  
  # get marginal likelihoods
  marglik1 <- list(
    logml = log_posterior(parameters=data.frame("r"=mean(priorr),"s"=mean(priors),"k"=mean(priork),"sigma.e"=mean(psige)
                                                ,"sigma.u"=mean(psigu)))
  )

  class(marglik1) <- "bridge"
  marglik1 <- JAGS_bridgesampling(fit=fit1, data=chimera, prior_list=priors_list, log_posterior)
  prob1 <- prior_prob(parameters=data.frame("r"=mean(priorr),"s"=mean(priors),"k"=mean(priork),"sigma.e"=mean(psige)
                                            ,"sigma.u"=mean(psigu)))
  
  bf <- inclusion_BF(prior_probs=c(prob0,prob1),margliks=c(as.numeric(marglik0[[1]]),as.numeric(marglik1[[1]])),is_null=c(T,F))
  estFit1$bf <- bf
  return(list(estFit1))                          
  
})
stopCluster(cl)   # stop the clusters
end <- Sys.time()
end - start









    #****************************model*******************

    # set up model
      mod = function(){
        #priors
        r~dgamma(1.25,5)    
        k~dgamma(11.25,0.75)   
        s~dgamma(1.8,0.6) 
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

  
  
    #********************************************************************
      runModel <- function(inits=list("r"=0.08,"k"=15,"s"=2.5)
                           ,params=c("tau.u","tau.e","r","k","s", "logTiter_pred")
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
    
    
    write.csv(postSamples,here(".\\inVitro\\JensUpdateLogisticFunction\\postDist1stSet_1.csv"))  
   write.csv(postChains,here(".\\inVitro\\JensUpdateLogisticFunction\\postDist1stSet_1_Ch.csv"))
  
   
    summaryStat <- lapply(unique(postSamples1$Chimera),function(X){

     temp2 <- postSamples2[postSamples2$Chimera %in% x,]
     r2Mean <- mean(temp2$r)
     s2Mean <- mean(temp2$s)
     k2Mean <- mean(temp2$k)

     return(c(pdr,pds,pdk))
   })
   
  

   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
  
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


  rPlot <- ggplot(data = postSamples, aes(x = r, y = Chimera)) + 
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


  kPlot <- ggplot(data = postSamples, aes(x = k, y = Chimera)) + 
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



  sPlot <- ggplot(data = postSamples, aes(x = s, y = Chimera)) + 
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

  sPlot
  rPlot
  kPlot



