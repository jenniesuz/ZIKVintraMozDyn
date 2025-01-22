library(ggridges)
library(ggmcmc)
library(grid)
library(gridExtra)
library(here)



#**************Read in the sensitivity results***********************
postSamples1Ch <- read.csv(here(".\\inVitro\\JensUpdateLogisticFunction\\postDist1stSet_1Ch.csv"))


#*************************chains:****************************
postSamples1Ch$prior <- "1"
postSamples2Ch$prior <- "2"
postSamples3Ch$prior <- "3"
postSamples4Ch$prior <- "4"
postSamples5Ch$prior <- "5"
postSamples6Ch$prior <- "6"
postSamples6Ch$prior <- "7"

chainData <- rbind.data.frame(postSamples1Ch,postSamples2Ch,postSamples3Ch,postSamples4Ch,postSamples5Ch,postSamples6Ch)
chainData$Chimera <- factor(chainData$Chimera
                            ,labels=c("1: Senegal"
                                      ,"2: Thai"
                                      ,"3: Senegal/Thai-SP"
                                      ,"4: Thai/Senegal-SP"
                                      ,"5: Senegal/Thai-nSP"
                                      ,"6: Thai/Senegal-nSP"
                                      ,"7: Senegal/Thai-UTR"
                                      ,"8: Thai/Senegal-UTR"))


# Use ggplot() to construct a trace plot 
pdf(file = here(".\\inVitro\\JensUpdateLogisticFunction\\figures\\rChains.pdf"),   # The directory you want to save the file in
    width = 20, # The width of the plot in inches
    height = 20) # The height of the plot in inches
ggplot(chainData, aes(x = iter, y = r)) + 
  geom_line() + theme_ridges() +
  facet_wrap(~ prior+Chimera)
dev.off()

pdf(file = here(".\\inVitro\\JensUpdateLogisticFunction\\figures\\kChains.pdf"))   # The directory you want to save the file in
ggplot(chainData, aes(x = iter, y = k)) + 
  geom_line() + theme_ridges() +
  facet_wrap(~prior+Chimera)
dev.off()

pdf(file = here(".\\inVitro\\JensUpdateLogisticFunction\\figures\\sChains.pdf"))   # The directory you want to save the file in
ggplot(chainData, aes(x = iter, y = s)) + 
  geom_line() + theme_ridges() +
  facet_wrap(~prior+Chimera)
dev.off()

rm(list=ls())
#****************posterior distributions******************
library(ggridges)
library(ggmcmc)
library(grid)
library(gridExtra)
library(here)

postSamples1 <- read.csv(here(".\\inVitro\\JensUpdateLogisticFunction\\postDist1stSet-r0.05_1-k5_0.5-s1.5_0.5.csv"))
postSamples2 <- read.csv(here(".\\inVitro\\JensUpdateLogisticFunction\\postDist1stSet-r0.05_1-k5_0.5-s3_1.csv"))
postSamples3 <- read.csv(here(".\\inVitro\\JensUpdateLogisticFunction\\postDist1stSet-r0.05_1-k10_1-s1.5_0.5.csv"))
postSamples4 <- read.csv(here(".\\inVitro\\JensUpdateLogisticFunction\\postDist1stSet-r0.05_1-k10_1-s3_1.csv"))
postSamples5 <- read.csv(here(".\\inVitro\\JensUpdateLogisticFunction\\postDist1stSet-r0.05_1-k10_1-s3_1.csv"))
postSamples6 <- read.csv(here(".\\inVitro\\JensUpdateLogisticFunction\\postDist1stSet-r1.25_50-k100_10-s9_3.csv"))
postSamples7 <- read.csv(here(".\\inVitro\\JensUpdateLogisticFunction\\postDist1stSet-r0.125_5-k100_10-s9_3.csv"))

postSamples1 <- postSamples1[,c("X","k","r","s","Chimera")]
postSamples3 <- postSamples3[,c("X","k","r","s","Chimera")]


postSamples1$prior <- "1"
postSamples2$prior <- "2"
postSamples3$prior <- "3"
postSamples4$prior <- "4"
postSamples5$prior <- "5"
postSamples6$prior <- "6"
postSamples7$prior <- "7"
post <- rbind.data.frame(postSamples1,postSamples2,postSamples3,postSamples4,postSamples5,postSamples6,postSamples7)

post$Chimera <- factor(post$Chimera
                       ,labels=c("1: Senegal"
                                 ,"2: Thai"
                                 ,"3: Senegal/Thai-SP"
                                 ,"4: Thai/Senegal-SP"
                                 ,"5: Senegal/Thai-nSP"
                                 ,"6: Thai/Senegal-nSP"
                                 ,"7: Senegal/Thai-UTR"
                                 ,"8: Thai/Senegal-UTR"))


pdf(file = here(".\\inVitro\\JensUpdateLogisticFunction\\figures\\kPosterior.pdf"))   # The directory you want to save the file in
ggplot(data = post, aes(x = k, y = Chimera)) + 
  stat_density_ridges(quantile_lines = TRUE, 
                      quantiles = c(0.025, 0.5, 0.975), 
                      alpha = 0.2
                      ,fill="blue") + 
  labs(title="Carrying capacity (k)") +
  xlab("Posterior estimate") + 
  ylab("") +
  theme_bw() +
  theme(panel.border = element_blank()                   
        ,axis.line = element_line(color = 'black')
        ,text=element_text(size=10)
        ,plot.margin=unit(c(0.2,0.1,0.1,0.1), "cm")
        ,axis.text=element_text(size=8)
        ,legend.key.size = unit(0.5,"line")  
        ,legend.background = element_blank()
        ,legend.text=element_text(size=8)
        ,strip.text=element_text(size=8)
        ,strip.background =element_blank()
  ) +
  facet_wrap(~prior)
dev.off()

pdf(file = here(".\\inVitro\\JensUpdateLogisticFunction\\figures\\sPosterior.pdf"))   # The directory you want to save the file in
ggplot(data = post, aes(x = s, y = Chimera)) + 
  stat_density_ridges(quantile_lines = TRUE, 
                      quantiles = c(0.025, 0.5, 0.975), 
                      alpha = 0.2
                      ,fill="blue") + 
  labs(title="Starting virus (s)") +
  xlab("Posterior estimate") + 
  ylab("") +
  theme_bw() +
  theme(panel.border = element_blank()                   
        ,axis.line = element_line(color = 'black')
        ,text=element_text(size=10)
        ,plot.margin=unit(c(0.2,0.1,0.1,0.1), "cm")
        ,axis.text=element_text(size=8)
        ,legend.key.size = unit(0.5,"line")  
        ,legend.background = element_blank()
        ,legend.text=element_text(size=8)
        ,strip.text=element_text(size=8)
        ,strip.background =element_blank()
  ) +
  facet_wrap(~prior)
dev.off()

pdf(file = here(".\\inVitro\\JensUpdateLogisticFunction\\figures\\rPosterior.pdf"))   # The directory you want to save the file in
ggplot(data = post, aes(x = r, y = Chimera)) + 
  stat_density_ridges(quantile_lines = TRUE, 
                      quantiles = c(0.025, 0.5, 0.975), 
                      alpha = 0.2
                      ,fill="blue") + 
  labs(title="Growth rate (r)") +
  xlab("Posterior estimate") + 
  ylab("") +
  theme_bw() +
  theme(panel.border = element_blank()                   
        ,axis.line = element_line(color = 'black')
        ,text=element_text(size=10)
        ,plot.margin=unit(c(0.2,0.1,0.1,0.1), "cm")
        ,axis.text=element_text(size=8)
        ,legend.key.size = unit(0.5,"line")  
        ,legend.background = element_blank()
        ,legend.text=element_text(size=8)
        ,strip.text=element_text(size=8)
        ,strip.background =element_blank()
  ) +
  facet_wrap(~prior)
dev.off()


#**************percentage deviation from 1*************

# compare results between 1,2 & 3 in particular - others seemed to have more problems with convergence
# by chimera
compareAv <- lapply(unique(postSamples1$Chimera),function(X){
  temp1 <- postSamples1[postSamples1$Chimera %in% x,]
  r1Mean <- mean(temp1$r)
  s1Mean <- mean(temp1$s)
  k1Mean <- mean(temp1$k)
  
  temp2 <- postSamples2[postSamples2$Chimera %in% x,]
  r2Mean <- mean(temp2$r)
  s2Mean <- mean(temp2$s)
  k2Mean <- mean(temp2$k)
  
  pdr <- (r1Mean-r2Mean)/r1Mean*100
  pds <- (s1Mean-s2Mean)/s1Mean*100
  pdk <- (k1Mean-k2Mean)/k1Mean*100
  return(c(pdr,pds,pdk))
})


lapply(unique(postSamples1$Chimera),function(x){
  temp<- postSamples1[postSamples1$Chimera %in% x,] 
mean(temp$s)})
