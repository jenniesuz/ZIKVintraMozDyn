


library(ggplot2)
library(gridExtra)
library(deSolve)
#*********************PARAMETERS*************************************
params <- function(           
    infRate = 0.005            # probability of contact and infection 
    ,scalingParameter = 0.5        
    ,cMax = 6       # density of cells
)
  return(as.list(environment()))
#***************************************************************

#*****************INITIAL CONDITIONS*****************************
initial <- c(V=1)

times <- seq(0,96,1)               # times to solve at
#**************************************************************

#****************MODEL*****************************************
mod <- function(tt,yy,parms) with(c(parms,as.list(yy)), {
  
  deriv <- infRate*V*(cMax-V)/(1+scalingParameter*(cMax-V)) 
  
  return(list(deriv))
})
#*************************************************************

#**************SIMULATE***************************************
simPop <- function(init=initial, tseq = times, modFunction=mod, parms = params()) {
  simDat <- as.data.frame(ode(init, tseq, modFunction, parms=parms))
  return(simDat)
}
#****************************************************************

library(here)
source(here(".//inVitro//JensUpdateExponentialFunction//dataSecondSet.R"))


c1 <- dat[dat$Chimera %in% "Chimera 9",]

sim <- simPop(parms=params(infRate=0.03
  ,scalingParameter=0.35
  ,cMax = 14
  ))

#plot(sim$time,sim$V)

ggplot(c1) +
  geom_point(aes(x=time,y=log(Titer))) +
  ylim(1,15) +
  geom_line(data=sim,aes(x=time,y=V))
