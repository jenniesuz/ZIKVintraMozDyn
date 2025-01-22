
library(ggplot2)
library(gridExtra)
library(deSolve)

#*********************PARAMETERS*************************************
params <- function(     
    muV =  0.2 #0.1          # virus clearance/ death rate 
    ,probInf = 10^-10        # probability of contact and infection 
    ,prodRate =   2          # virion production rate
    ,cellSpread = 0          # rate virions in infected cells infect susceptible cells
    ,cMax = 10^6             # number of cells in midgut
)
  return(as.list(environment()))
#***************************************************************

#*****************INITIAL CONDITIONS*****************************
initial <- c(Gv = 10^5      # number of virions in bloodmeal
             ,Mci = 0       # number of infected midgut cells
             ,Mv = 0        # number of virions in midgut
)

times <- seq(0,96,1)               # times to solve at
#**************************************************************

#****************MODEL*****************************************
mod <- function(tt,yy,parms) with(c(parms,as.list(yy)), {
  
  deriv <- rep(NA,3)
  
  deriv[1] <- -Gv*probInf*(cMax - Mci) - muV*Gv      # virion in bloodmeal
  
  deriv[2] <- Gv*probInf*(cMax-Mci) + cellSpread*Mci*(cMax-Mci)   # infected midgut cells
  
  deriv[3] <- prodRate*Mci - muV*Mv       # virions in midgut
  
  return(list(deriv))
})
#*************************************************************

#**************SIMULATE***************************************
simPop <- function(init=initial, tseq = times, modFunction=mod, parms = params()) {
  simDat <- as.data.frame(ode(init, tseq, modFunction, parms=parms))
  return(simDat)
}
#****************************************************************

sim <- simPop(parms=params())

ggplot(sim) +
  geom_line(aes(x=time,y=log10(Mv))) 

