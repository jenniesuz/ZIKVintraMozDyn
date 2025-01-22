
library(here)
source(here(".//inVitro//functionInfectionModelWithEclipseWithCellDeath.R"))
source(here(".//inVitro//functionLikelihood.R"))
source(here(".//inVitro//data.R"))
library(bbmle)


dataPlot +
  #geom_line(data=test,aes(x=Time,y=log10(PFU)),linetype=2) +
  geom_line(data=sim,aes(x=Time,y=log10(PFU)),linetype=3)


#* first try with one chimera
senegal <- datLong[(datLong$Chimera %in% "Senegal strain (Chimera 1)") & (datLong$Time %in% 24 | datLong$Time %in% 48),]

senegalFit <- mle2(function(par1,par2){nllPoisTwoPar(par1,par2,senegal)}                       
              ,start=list(par1=log(10^-8),par2=log(10)))

senegalCoef <- exp(coef(senegalFit))
senegalSim <- simPop(parms=params(probInf=senegalCoef[1],prodRate=senegalCoef[2]))


# separate fits for each - easiest approach to start with
# but will probably want to compare models assuming same/ different parameters

#*****************fits*****************
fits <- lapply(unique(datLong$Chimera),function(x){
  temp <- datLong[(datLong$Chimera %in% x) & (datLong$Time %in% 24 | datLong$Time %in% 48),]
  
  fit <- mle2(function(par1,par2){nllPois(par1,par2,temp)}                       
                     ,start=list(par1=log(10^-8),par2=log(10)))
  
  coefs <- exp(coef(fit))
})

fits <- do.call(rbind.data.frame,fits)
fits$Chimera <- unique(datLong$Chimera)
names(fits) <- c("probInf","prodRate","Chimera")


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
