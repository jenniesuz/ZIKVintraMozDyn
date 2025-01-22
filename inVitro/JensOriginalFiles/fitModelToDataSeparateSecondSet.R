library(here)
source(here(".//inVitro//JensOriginalFiles//functionInfectionModelWithEclipse.R"))
source(here(".//inVitro//JensOriginalFiles//functionLikelihood.R"))
source(here(".//inVitro//JensUpdateUsingMikesCode//dataSecondSet.R"))
library(bbmle)

#* first try with one chimera
c1 <- dat[dat$Chimera %in% "Chimera 1",]


simTest <- simPop(parms=params(muV=2
                               ,probInf=10^-2
                               ,prodRate=20
                               ,cellSpread=0
                               ,cMax=10^6
                               ,eclipse=1/24))

ggplot(c1) +
  geom_point(aes(x=time
                 ,y=log10(Titer)
  )) +
  geom_line(data=simTest,aes(x=time,y=log10(Mv))) +
  theme_grey()






c1Fit <- mle2(function(par1,par2){nllPois(par1,par2,senegal)}                       
              ,start=list(par1=log(10^-8),par2=log(10)))

senegalCoef <- exp(coef(senegalFit))
senegalSim <- simPop(parms=params(probInf=senegalCoef[1],prodRate=senegalCoef[2]))


# separate fits for each - easiest approach to start with
# but will probably want to compare models assuming same/ different parameters

#*****************fits*****************
params <- function(     
    muV =  0.001 #0.1          # virus clearance/ death rate 
    ,probInf = 10^-8           # probability of contact and infection 
    ,prodRate = 15             # virion production rate
    ,cMax = 10^6               # number of cells in midgut
    ,eclipse = 1/24
    ,cDeath =1/36
)
  return(as.list(environment()))


#*******Only inf rate fitted********
infRateFits <- lapply(unique(datLong$Chimera),function(x){
  temp <- datLong[(datLong$Chimera %in% x) & (datLong$Time %in% 24 | datLong$Time %in% 48),]
  
  fit <- mle2(function(par1){nllPoisInfRate(par1,temp)}                       
              ,start=list(par1=log(10^-8)))
  
  coefs <- exp(coef(fit))
  return(c(coefs,AIC(fit)))
})

infRateFits <- do.call(rbind.data.frame,infRateFits)
infRateFits$Chimera <- unique(datLong$Chimera)
names(infRateFits) <- c("probInf","AIC","Chimera")


#********prod rate fits********
prodRateFits <- lapply(unique(datLong$Chimera),function(x){
  temp <- datLong[(datLong$Chimera %in% x) & (datLong$Time %in% 24 | datLong$Time %in% 48),]
  
  fit <- mle2(function(par1){nllPoisProdRate(par1,temp)}                       
              ,start=list(par1=log(20)))
  
  coefs <- exp(coef(fit))
  return(c(coefs,AIC(fit)))
})

prodRateFits <- do.call(rbind.data.frame,prodRateFits)
prodRateFits$Chimera <- unique(datLong$Chimera)
names(prodRateFits) <- c("probInf","AIC","Chimera")



#*************Both parameters fitted****************
fits <- lapply(unique(datLong$Chimera),function(x){
  temp <- datLong[(datLong$Chimera %in% x) & (datLong$Time %in% 24 | datLong$Time %in% 48),]
  
  fit <- mle2(function(par1,par2){nllPoisTwoPar(par1,par2,temp)}                       
                     ,start=list(par1=log(10^-8),par2=log(10)))
  
  coefs <- exp(coef(fit))
  return(c(coefs,AIC(fit)))
})

fits <- do.call(rbind.data.frame,fits)
fits$Chimera <- unique(datLong$Chimera)
names(fits) <- c("probInf","prodRate","AIC","Chimera")






#*************sims*************
sims <- lapply(unique(fits$Chimera),function(x){
  temp <- fits[fits$Chimera %in% x,]
   sim <- simPop(parms=params(probInf=temp[1],prodRate=temp[2]))
   sim$Chimera <- x
return(sim)
})
sims <- do.call(rbind.data.frame,sims)
names(sims)[4] <- "PFU"



#*********************PLOT*******************

datLong <- datLong[!datLong$Time %in% 72,]
datLong <- datLong[!datLong$Time %in% 96,]
sims <- sims[!sims$Time >50,]

ggplot(datLong) +
  geom_point(aes(x=Time 
                ,y=PFU
                ,group=Replicate)) +
  geom_line(data=sims,aes(x=Time,y=PFU),linetype=2) +
  xlim(24,50) +
 # ylim(0,10^6) +
  scale_color_manual("Chimera", values = cols) +
  theme_grey() +
  facet_wrap(~Chimera)
