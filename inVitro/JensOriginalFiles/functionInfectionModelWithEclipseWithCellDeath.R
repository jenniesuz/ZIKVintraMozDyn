
library(ggplot2)
library(gridExtra)
library(deSolve)

#*********************PARAMETERS*************************************
params <- function(     
    muV =  0.001 #0.1          # virus clearance/ death rate 
    ,probInf = 10^-8           # probability of contact and infection 
    ,prodRate = 15             # virion production rate
    ,cMax = 10^6               # number of cells in midgut
    ,eclipse = 1/24
    ,cDeath =1/36
)
  return(as.list(environment()))
#***************************************************************

#*****************INITIAL CONDITIONS*****************************
initial <- c(Me = 5
             ,Mci = 0       # number of infected midgut cells
             ,Mv = 0        # number of virions in midgut
)

times <- seq(0,100,1)               # times to solve at
#**************************************************************

#****************MODEL*****************************************
mod <- function(tt,yy,parms) with(c(parms,as.list(yy)), {
  
  deriv <- rep(NA,3)
  
  deriv[1] <- Mv*probInf*(cMax - Mci - Me) - eclipse*Me - Me*cDeath 
  
  deriv[2] <- eclipse*Me  - Mci*cDeath # infected midgut cells
  
  deriv[3] <- prodRate*Mci - muV*Mv       # virions in midgut
  
  return(list(deriv))
})
#*************************************************************

#**************SIMULATE***************************************
simPop <- function(init=initial, tseq = times, modFunction=mod, parms = params()) {
  simDat <- as.data.frame(ode(init, tseq, modFunction, parms=parms))
  names(simDat)[1] <- "Time"
  return(simDat)
}
#****************************************************************

sim <- simPop(parms=params())

ggplot(sim) +
  geom_line(aes(x=Time,y=log10(Mv)))

names(sim)[4]<-"PFU"

