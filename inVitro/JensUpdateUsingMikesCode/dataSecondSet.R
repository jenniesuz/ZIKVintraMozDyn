library(here)
library(ggplot2)

dat <- read.csv(here(".//data//growth2ndChimeras0.01JSL.csv"))

head(dat)


dataPlot <- ggplot(dat) +
  geom_point(aes(x=time
                ,y=log10(Titer)
                )) +
  facet_wrap(~Chimera) +
  xlim(0,100) +
  ylim(0,6) +
  theme_grey()

dataPlot


ggplot(dat) +
  geom_point(aes(x=time
                 ,y=Titer
  )) +
  facet_wrap(~Chimera) +
  xlim(0,100) +
  #ylim(0,6) +
  theme_grey()
