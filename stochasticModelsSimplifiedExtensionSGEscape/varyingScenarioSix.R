library(ggplot2)
#**************88gamma dists*********
oneMean = 0.00008
oneVar = 10^-7.5

oneSc <- oneVar/oneMean
oneSh <- oneMean^2/oneVar

oneSh <- 0.9
oneSc <- 0.0001
x.max = qgamma(0.999, shape=oneSh, scale=oneSc)
one.x = seq(from=0, to=x.max, by=x.max/1000)
one.dens = dgamma(one.x, shape=oneSh, scale=oneSc)
plot(one.x, one.dens, type='l')

rgamma(10, shape=oneSh, scale = oneSc)

twoMean = 0.00005
twoVar = 10^-7.3

twoSc <- twoVar/twoMean
twoSh <- twoMean^2/twoVar
rate <- 1/twoSc
x.max = qgamma(0.999, shape=twoSh, scale=twoSc)
two.x = seq(from=0, to=x.max, by=x.max/1000)
two.dens = dgamma(two.x, shape=twoSh, scale=twoSc)
plot(two.x, two.dens, type='l')



threeMean = 0.00005
threeVar = 10^-7

threeSc <- threeVar/threeMean
threeSh <- threeMean^2/threeVar
rate <- 1/threeSc
x.max = qgamma(0.999, shape=threeSh, scale=threeSc)
three.x = seq(from=0, to=x.max, by=x.max/1000)
three.dens = dgamma(three.x, shape=threeSh, scale=threeSc)
plot(three.x, three.dens, type='l')

fourMean = 0.00003
fourVar = 10^-6.5

fourSc <- fourVar/fourMean
fourSh <- fourMean^2/fourVar
rate <- 1/fourSc
x.max = qgamma(0.999, shape=fourSh, scale=fourSc)
four.x = seq(from=0, to=x.max, by=x.max/1000)
four.dens = dgamma(four.x, shape=fourSh, scale=fourSc)
plot(four.x, four.dens, type='l')

dist <- c(rep("1",length(one.x))
          ,rep("2",length(two.x))
          ,rep("3",length(three.x))
          ,rep("4",length(four.x)))

xvals <- c(one.x,two.x,three.x,four.x)
yvals <- c(one.dens,two.dens,three.dens,four.dens)

dists <- cbind.data.frame(dist,xvals,yvals)

ggplot(dists) +
  geom_line(aes(x=log10(xvals),y=yvals,col=dist,group=dist)) 

#*************varing variation bewteen mosquitoes in h:s escape********

modelOutput <- repeatModel(virus_params(bloodmealClearance = 1/72
                                        ,propSuccessInf = 10^-4
                                        ,growthRateH = 0.04
                                        ,growthRateM = 0.04
                                        ,carryCap = 10^20
                                        ,escapeRateM =  0.00005 
                                        ,escapeRateH = 0.00008
                                        ,escapeVarS = 10^-7.5
                                        ,escapeVarSTrue = T) # 0.001, 0.00005
                           ,startingVirus=10^6)

modelOutput <- do.call(rbind,modelOutput)

propsSen8a <- dissSummaryFunc(modelOutput)
propsSen8a$run <- "1"

#***************
modelOutput <- repeatModel(virus_params(bloodmealClearance = 1/72
                                        ,propSuccessInf = 10^-4
                                        ,growthRateH = 0.04
                                        ,growthRateM = 0.04
                                        ,carryCap = 10^20
                                        ,escapeRateM = 0.00005 
                                        ,escapeRateH = 0.00005
                                        ,escapeVarS = 10^-7.3
                                        ,escapeVarSTrue = T) # 0.001, 0.00005
                           ,startingVirus=10^6)

modelOutput <- do.call(rbind,modelOutput)

propsSen8b <- dissSummaryFunc(modelOutput)
propsSen8b$run <- "2"


#**************************
modelOutput <- repeatModel(virus_params(bloodmealClearance = 1/72
                                        ,propSuccessInf = 10^-4
                                        ,growthRateH = 0.04
                                        ,growthRateM = 0.04
                                        ,carryCap = 10^20
                                        ,escapeRateM = 0.00005 
                                        ,escapeRateH = 0.00005
                                        ,escapeVarS = 10^-7
                                        ,escapeVarSTrue = T) # 0.001, 0.00005
                           ,startingVirus=10^6)

modelOutput <- do.call(rbind,modelOutput)

propsSen8c <- dissSummaryFunc(modelOutput)
propsSen8c$run <- "3"



modelOutput <- repeatModel(virus_params(bloodmealClearance = 1/72
                                        ,propSuccessInf = 10^-4
                                        ,growthRateH = 0.04
                                        ,growthRateM = 0.04
                                        ,carryCap = 10^20
                                        ,escapeRateM = 0.00005 
                                        ,escapeRateH = 0.00003
                                        ,escapeVarS = 10^-6.5
                                        ,escapeVarSTrue = T) # 0.001, 0.00005
                           ,startingVirus=10^6)

modelOutput <- do.call(rbind,modelOutput)

propsSen8d <- dissSummaryFunc(modelOutput)
propsSen8d$run <- "4"



propsSen <- rbind.data.frame(propsSen8a,propsSen8b,propsSen8c,propsSen8d)
names(propsSen)[9] <- "Virus"
#******************
#*
#*
cols <- c("#8c510a"
  ,"#bf812d"
  ,"#dfc27d"
  ,"#f6e8c3"
  ,"#c7eae5"
  ,"#80cdc1"
  ,"#35978f"
  ,"#01665e"
)

cols <- c("#8c510a"
          ,"#dfc27d"
          ,"#c7eae5"
          ,"#01665e")





tiff(here::here(".//stochasticModelsSimplifiedExtensionSGEscape//fig_simsChimerasTrans.tiff")
     , height =6, width = 4, units = 'in', compression="lzw", res=400)

ggplot(propsSen) +
  #geom_point(data=dat2,aes(x=time,y=propTrans,col=Virus),alpha=0.5) +
  #geom_errorbar(data=dat2,aes(x=time,ymin=propTransLower,ymax=propTransupper,col=Virus),alpha=0.5) +
  geom_line(aes(x=time,y=proportionCapableofTransmission,col=Virus,group=Virus),linetype="longdash") +
  scale_color_manual(name="Virus",values=cols) +
  xlim(5,16) +
  xlab("Time (days)") +
  labs(col="Variance (log10)") +
  ylab("Proportion with virus \n in salivary glands") +
  theme_set(theme_bw())  +
  theme(panel.border = element_blank()
        ,axis.line = element_line(color = 'black')
        ,text=element_text(size=12)
        ,plot.margin=unit(c(0.2,0.3,0.1,0.1), "cm")
        ,axis.text=element_text(size=10)
        ,legend.key.size = unit(0.5,"line")
        ,legend.background = element_blank()
        ,legend.text=element_text(size=10)
        ,legend.position ="none"
        ,legend.title = element_blank()
  )

dev.off()





tiff(here::here(".//stochasticModelsSimplifiedExtensionSGEscape//fig_simsChimerasInf.tiff")
     , height =6, width = 4, units = 'in', compression="lzw", res=400)

ggplot(propsSen) +
  #geom_point(data=dat2,aes(x=time,y=propInf,col=as.factor(Virus)),alpha=0.5) +
  #geom_errorbar(data=dat2,aes(x=time,ymin=propInfLower,ymax=propInfupper,col=Virus),alpha=0.5) +
  geom_line(aes(x=time,y=proportionInfected,col=Virus,group=Virus),linetype="longdash") +
  scale_color_manual(name="Run",values=cols) +
  xlim(5,16) +
  xlab("Time (days)") +
  ylab("Proportion infected") +
  labs(col="Run") +
  theme_set(theme_bw())  +
  theme(panel.border = element_blank()
        ,axis.line = element_line(color = 'black')
        ,text=element_text(size=12)
        ,plot.margin=unit(c(0.2,0.3,0.1,0.1), "cm")
        ,axis.text=element_text(size=10)
        ,legend.key.size = unit(0.5,"line")
        ,legend.background = element_blank()
        ,legend.text=element_text(size=10)
        ,legend.position =c(0.8,0.2)
        
  )
dev.off()







tiff(here::here(".//stochasticModelsSimplifiedExtensionSGEscape//fig_simsChimerasDiss.tiff")
     , height =6, width = 4, units = 'in', compression="lzw", res=400)

ggplot(propsSen) +
 # geom_point(data=dat2,aes(x=time,y=propDiss,col=Virus),alpha=0.5) +
#  geom_errorbar(data=dat2,aes(x=time,ymin=propDissLower,ymax=propDissupper,col=Virus),alpha=0.5) +
  geom_line(aes(x=time,y=proportionDisseminated,col=Virus,group=Virus),linetype="longdash") +
  scale_color_manual(name="Virus",values=cols) +
  xlim(5,16) +
  xlab("Time (days)") +
  labs(col="Variance (log10)") +
  ylab("Proportion disseminated") +
  theme_set(theme_bw())  +
  theme(panel.border = element_blank()
        ,axis.line = element_line(color = 'black')
        ,text=element_text(size=12)
        ,plot.margin=unit(c(0.2,0.3,0.1,0.1), "cm")
        ,axis.text=element_text(size=10)
        ,legend.key.size = unit(0.5,"line")
        ,legend.background = element_blank()
        ,legend.text=element_text(size=10)
        ,legend.position ="none"
        ,legend.title = element_blank()
  )
dev.off()







saveRDS(propsSen,"varyingGammaDist.rds")


#***************observed data***********

dat <- read.csv(here(".\\data\\030823_Chimera_1-8 time to dissemination.csv"))
head(dat)

dat2 <- lapply(unique(dat$Virus),function(x){
  temp <- dat[dat$Virus %in% x,]
  seven <- temp[temp$Day %in% 7,]
  totSeven <- colSums(seven[5:8])
  sampSeven <- length(seven$Infection)
  
  ten <- temp[temp$Day %in% 10,]
  totTen<- colSums(ten[5:8])
  sampTen <- length(ten$Infection)
  
  fourt <- temp[temp$Day %in% 14,]
  totFourt<- colSums(fourt[5:8])
  sampFourt <- length(fourt$Infection)
  
  Virus <- rep(x,3)
  Samp <- c(sampSeven,sampTen,sampFourt)
  Infection <- c(totSeven[1],totTen[1],totFourt[1])
  Dissemination <- c(totSeven[2],totTen[2],totFourt[2])
  Transmission <- c(totSeven[3],totTen[3],totFourt[3])
  day <- c(7,10,14)
  return(cbind.data.frame(Virus,Samp,Infection,Dissemination,Transmission,day))
  
})

library(binom)
dat2 <- do.call(rbind,dat2)
inf <- binom.confint(x=dat2$Infection,n=dat2$Samp,method="exact")
diss <- binom.confint(x=dat2$Dissemination,n=dat2$Samp,method="exact")
trans <- binom.confint(x=dat2$Transmission,n=dat2$Samp,method="exact")

dat2$propInf <- inf$mean
dat2$propInfLower <- inf$lower
dat2$propInfupper <- inf$upper

dat2$propDiss <- diss$mean
dat2$propDissLower <- diss$lower
dat2$propDissupper <-diss$upper


dat2$propTrans <- trans$mean
dat2$propTransLower <- trans$lower
dat2$propTransupper <- trans$upper


names(dat2)[6] <-"time"

ggplot(dat2) +
  geom_point(aes(x=day,y=propInf))

#***


