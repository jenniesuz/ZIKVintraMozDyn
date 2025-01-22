library(here)
library(ggplot2)

dat <- read.csv(here(".//data//growth1stChimeras0.01ST240424.csv"))

head(dat)

dat$Chimera <- factor(dat$Chimera
                              ,labels=c("Senegal"
                                        ,"Thai"
                                        ,"Senegal/Thai-SP"
                                        ,"Thai/Senegal-SP"
                                        ,"Senegal/Thai-nSP"
                                        ,"Thai/Senegal-nSP"
                                        ,"Senegal/Thai-UTR"
                                        ,"Thai/Senegal-UTR"))


ggplot(dat) +
  geom_point(aes(x=time
                 ,y=log10(Titer)
  )) +
  #xlim(0,100) +
  #ylim(0,20) +
  facet_wrap(~Chimera) +
  theme_grey()

ggplot(dat) +
  geom_point(aes(x=time
                ,y=log(Titer)
                )) +
  xlim(0,100) +
  ylim(0,20) +
  facet_wrap(~Chimera) +
  theme_grey()





