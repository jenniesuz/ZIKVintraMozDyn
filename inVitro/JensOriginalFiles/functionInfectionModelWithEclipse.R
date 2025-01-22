
library(ggplot2)
library(gridExtra)
library(deSolve)

#*********************PARAMETERS*************************************
params <- function(     
    muV =  0.001 #0.1          # virus clearance/ death rate 
    ,probInf = 10^-12        # probability of contact and infection 
    ,prodRate =   5          # virion production rate
    ,cellSpread = 0          # rate virions in infected cells infect susceptible cells
    ,cMax = 10^6             # number of cells in midgut
    ,eclipse = 1/24
)
  return(as.list(environment()))
#***************************************************************

#*****************INITIAL CONDITIONS*****************************
initial <- c(Gv = 0     # number of virions in bloodmeal
             ,Me = 10^6*0.01
             ,Mci = 0       # number of infected midgut cells
             ,Mv = 0        # number of virions in midgut
)

times <- seq(0,100,1)               # times to solve at
#**************************************************************

#****************MODEL*****************************************
mod <- function(tt,yy,parms) with(c(parms,as.list(yy)), {
  
  deriv <- rep(NA,4)
  
  deriv[1] <- -Gv*probInf*(cMax - Mci - Me) - muV*Gv      # virion in bloodmeal
  
  deriv[2] <- Gv*probInf*(cMax - Mci - Me) + cellSpread*Mci*(cMax-Mci-Me) - eclipse*Me 
  
  deriv[3] <- eclipse*Me    # infected midgut cells
  
  deriv[4] <- prodRate*Mci - muV*Mv       # virions in midgut
  
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
  geom_line(aes(x=time,y=log10(Mv))) # +
 # xlim(24,50)

