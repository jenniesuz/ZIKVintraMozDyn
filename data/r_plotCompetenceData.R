library(here)
library(ggplot2)


dat <- read.csv(here(".//data//30052023_Chimera_1-8_dose_response_curve_JSL.csv"))


#************Dose response curve*************
library(plyr)
library(binom)

infPrev <- ddply(dat,.(Chimera,Dose),summarise,denom=length(ID),num=sum(RT.PCR))
infPrevb <- binom.confint(x=infPrev$num,n=infPrev$denom,methods="exact")
infPrev$mean <- infPrevb$mean
infPrev$lower <- infPrevb$lower
infPrev$upper <- infPrevb$upper

infPrev$Chimera <- factor(infPrev$Chimera
                          ,labels=c("1: Senegal"
                                    ,"2: Thai"
                                    ,"3: Senegal/Thai-SP"
                                    ,"4: Thai/Senegal-SP"
                                    ,"5: Senegal/Thai-nSP"
                                    ,"6: Thai/Senegal-nSP"
                                    ,"7: Senegal/Thai-UTR"
                                    ,"8: Thai/Senegal-UTR"))



infPrev$ChimeraNSPTN <- infPrev$Chimera  
infPrev$ChimeraNSPTN <- factor(infPrev$ChimeraNSPTN
                               ,labels=c("Senegal-nSP"
                                         ,"Thai-nSP"
                                         ,"Senegal-nSP"
                                         ,"Thai-nSP"
                                         ,"Thai-nSP"
                                         ,"Senegal-nSP"
                                         ,"Senegal-nSP"
                                         ,"Thai-nSP"))



infPrev$ChimeraSPTN <- infPrev$Chimera  
infPrev$ChimeraSPTN <- factor(infPrev$ChimeraSPTN
                              ,labels=c("Senegal-SP"
                                        ,"Thai-SP"
                                        ,"Thai-SP"
                                        ,"Senegal-SP"
                                        ,"Senegal-SP"
                                        ,"Thai-SP"
                                        ,"Senegal-SP"
                                        ,"Thai-SP"))




ggplot(infPrev) +
  geom_point(aes(x=log10(Dose),y=mean,col=Chimera)) +
  geom_errorbar(aes(x=log10(Dose),ymin=lower,ymax=upper,col=Chimera)) +
  xlim(0,10) +
  ylab("Proportion of mosquitoes infected")

ggplot(infPrev) +
  geom_point(aes(x=log10(Dose),y=mean,col=ChimeraNSPTN)) +
  geom_errorbar(aes(x=log10(Dose),ymin=lower,ymax=upper,col=ChimeraNSPTN)) +
  xlim(2.5,7.5) +
  ylab("Proportion of mosquitoes infected") 



ggplot(infPrev) +
  geom_point(aes(x=log10(Dose),y=mean,col=ChimeraSPTN)) +
  geom_errorbar(aes(x=log10(Dose),ymin=lower,ymax=upper,col=ChimeraSPTN)) +
  xlim(2.5,7.5) +
  ylab("Proportion of mosquitoes infected") 