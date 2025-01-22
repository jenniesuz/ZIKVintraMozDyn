
library(deSolve)

vir.kin.mod<-function(t,y,params){ 
  U <- y[1] # uninfected target host cells
  E <- y[2]
  I <- y[3] # infected productively infected cells
  V <- y[4] # infectious viral titer
  
  # parameters
  beta<-params[1] # rate constant characterizing infection
  p<-params[2] # rate of viral titer increase per cell
  c<-params[3] # viral clearance rate
  e<-params[4] # eclipse
  mu<-params[5] # cell death
  
  # ODE's
  dU.dt <- -beta*U*V - U*mu
  dE.dt <- beta*U*V - e*E - E*mu
  dI.dt <- e*E - I*mu
  dV.dt <- p*I - c*V
  # list containing the derivatives
  xdot <- c(dU.dt,dE.dt,dI.dt,dV.dt)
  return(list(xdot))
}


vir.kin.sim<-function(beta,p,c,e,mu,v0,u0,t.max=100,ts=1){
  U0<-u0
  E0<-0
  I0<-0
  V0<-v0
  # timesteps
  times<-seq(0,t.max,ts) 
  # initial conditions
  init<-c(uninf=U0,exp=E0,inf=I0,vir=V0)
  return(data.frame(lsoda(init,times,vir.kin.mod,c(beta,p,c,e,mu)))) 
}

v <- 10^3 
u <- 10^6
inf.cons <- 10^-8 #10^-10# beta
vir.rate <- 10 # p
clear.rate <- 0.001 # c
eclipse <- 1/12
cellDeath <- 1/36

test<-vir.kin.sim(beta=inf.cons
                  ,p=vir.rate
                  ,c=clear.rate
                  ,e=eclipse
                  ,mu=cellDeath
                  ,v0=v
                  ,u0=u
                  ,t.max=100
                  ,ts=1) 

names(test)[5] <- "PFU"
names(test)[1] <- "Time"

dataPlot +
  geom_line(data=test,aes(x=Time,y=log10(PFU)),linetype=2)

#ggplot(test) +
#  geom_line(aes(x=Time,y=log10(PFU)),linetype=2)

