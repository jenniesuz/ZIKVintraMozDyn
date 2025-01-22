library(here)
source(here(".//inVitro//JensUpdateUsingMikesCode//model.R"))
source(here(".//inVitro//JensUpdateUsingMikesCode//likelihoodFunction.R"))
source(here(".//inVitro//JensUpdateUsingMikesCode//dataSecondSet.R"))
library(bbmle)

#* first try with one chimera
c1 <- dat[dat$Chimera %in% "Chimera 1",]

simTest <- simPop(parms=params(infRate=0.045
                               ,muV=0.02
                               ,scalingParameter=0.5))


ggplot(c1) +
  geom_point(aes(x=time
                 ,y=Titer
  )) +
  geom_line(data=simTest,aes(x=time,y=V)) +
  theme_grey() 




c1Fit <- mle2(function(par1,par2,par3,par4){nll(par1,par2,par3,par4,c1)}                       
              ,start=list(par1=log(0.046),par2=log(0.02),par3=log(0.5),par4=log(2)))

c1Coef <- exp(coef(c1Fit))
c1Coef[1]
c1Sim <- simPop(parms=params(infRate=c1Coef[1]
                             ,muV=c1Coef[2]
                             ,scalingParameter=c1Coef[3]))


ggplot(c1) +
  geom_point(aes(x=time
                 ,y=Titer
  )) +
  geom_line(data=c1Sim,aes(x=time,y=V)) +
  theme_grey() 


## Without fitting mu
c1Fit <- mle2(function(par1,par3,par4){nll2(par1,par3,par4,c1)}                       
              ,start=list(par1=log(0.046),par3=log(0.5),par4=log(2)))

c1Coef <- exp(coef(c1Fit))
c1Coef[1]
c1Sim <- simPop(parms=params(infRate=c1Coef[1]
                             ,scalingParameter=c1Coef[2]))


ggplot(c1) +
  geom_point(aes(x=time
                 ,y=Titer
  )) +
  geom_line(data=c1Sim,aes(x=time,y=V)) +
  theme_grey() 



## 

fitFunc <- function(ChimeraName){

  c2 <- dat[dat$Chimera %in% ChimeraName,]


  c2Fit <- mle2(function(par1,par3,par4){nll2(par1,par3,par4,c2)}                       
              ,start=list(par1=log(0.046),par3=log(0.5),par4=log(2)))

  c2Coef <- exp(coef(c2Fit))
  c2Coef[1]
  c2Sim <- simPop(parms=params(infRate=c2Coef[1]
                             ,scalingParameter=c2Coef[2]))


  p <- ggplot(c2) +
    geom_point(aes(x=time
                   ,y=Titer
    )) +
    geom_line(data=c2Sim,aes(x=time,y=V)) +
    ylim(0, 1850000)
    theme_grey() 
  
  return(list(c2Coef,p))
  }



fitFunc(unique(dat$Chimera)[1])
fitFunc(unique(dat$Chimera)[2])
fitFunc(unique(dat$Chimera)[3])
fitFunc(unique(dat$Chimera)[4])
fitFunc(unique(dat$Chimera)[5])
fitFunc(unique(dat$Chimera)[6])
fitFunc(unique(dat$Chimera)[7])
fitFunc(unique(dat$Chimera)[8])
fitFunc(unique(dat$Chimera)[9])
