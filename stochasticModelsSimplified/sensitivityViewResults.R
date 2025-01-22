library(here)
library(ggplot2 )
library(reshape2)


#********************************************Midgut results**************************************************
midRes <- readRDS(here(".//stochasticModelsSimplified//sensitivityMidgut240917.rds"))


midRes <- melt(midRes,id.vars = c("proportionInfected","run"))



midRes$variable <- factor(midRes$variable
                        , levels = c("bmClear"
                                     ,"probabilityInfection"
                                     ,"growthRate"
                                     ,"carryCap"
                        )
                        , ordered = TRUE
                        , labels=c("Bloodmeal clearance rate"
                                   ,"Probability of infection"
                                   ,"Growth rate"
                                   ,"Carrying capacity")
                        
                                   )


tiff(here(".//stochasticModelsSimplified//fig_midgutSensitivity.tiff"), height =4, width = 5, units = 'in', compression="lzw", res=400)
ggplot(midRes) +
  geom_point(aes(x=value,y=proportionInfected),col="grey",alpha=0.5,size=0.5) +
  facet_wrap(~variable, scales="free") +
  xlab("Parameter value") +
  ylab("Proportion of simulations with midgut infection") +
  theme_bw() +
  theme(panel.border = element_blank()                   
        ,axis.line = element_line(color = 'black')
        ,text=element_text(size=10)
        ,plot.margin=unit(c(0.2,0.1,0.1,0.1), "cm")
        ,axis.text=element_text(size=6)
        ,legend.key.size = unit(0.5,"line")  
        ,legend.background = element_blank()
        ,legend.text=element_text(size=8)
        ,strip.text=element_text(size=8)
        ,strip.background =element_blank()
  )
dev.off()




#*******************************Hemocoel and salivary gland results**************************************
hcRes <- readRDS(here(".//stochasticModelsSimplified//sensitivityHS240917.rds"))
names(hcRes)[1:6] <- c("sevenDaysH","tenDaysH","fourteenDaysH","sevenDaysS","tenDaysS","fourteenDaysS")

hcRes <- melt(hcRes,id.vars = c("sevenDaysH","tenDaysH","fourteenDaysH","sevenDaysS","tenDaysS","fourteenDaysS","run"))

hcRes$variable <- factor(hcRes$variable
                          , levels = c("bmClear"
                                       ,"probabilityInfection"
                                       ,"growthRate"
                                       ,"carryCap"
                                       ,"escapeRate"
                          )
                          , ordered = TRUE
                          , labels=c("Bloodmeal clearance rate"
                                     ,"Probability of infection"
                                     ,"Growth rate"
                                     ,"Carrying capacity"
                                     ,"Escape rate")
                          
)


tiff(here(".//stochasticModelsSimplified//fig_HSevenDaysSensitivity2.tiff"), height =4, width = 5, units = 'in', compression="lzw", res=400)
ggplot(hcRes) +
  geom_point(aes(x=value,y=sevenDaysH),col="grey",alpha=0.5,size=0.5) +
  facet_wrap(~variable, scales="free") +
  xlab("Parameter value") +
  ylab("Proportion of simulations with dissemination") +
  theme_bw() +
  theme(panel.border = element_blank()                   
        ,axis.line = element_line(color = 'black')
        ,text=element_text(size=10)
        ,plot.margin=unit(c(0.2,0.1,0.1,0.1), "cm")
        ,axis.text=element_text(size=6)
        ,legend.key.size = unit(0.5,"line")  
        ,legend.background = element_blank()
        ,legend.text=element_text(size=8)
        ,strip.text=element_text(size=8)
        ,strip.background =element_blank()
  )
dev.off()


tiff(here(".//stochasticModelsSimplified//fig_SsevenDaysSensitivity2.tiff"), height =4, width = 5, units = 'in', compression="lzw", res=400)
ggplot(hcRes) +
  geom_point(aes(x=value,y=sevenDaysS),col="grey",alpha=0.5,size=0.5) +
  facet_wrap(~variable, scales="free") +
  xlab("Parameter value") +
  ylab("Proportion of simulations with dissemination") +
  theme_bw() +
  theme(panel.border = element_blank()                   
        ,axis.line = element_line(color = 'black')
        ,text=element_text(size=10)
        ,plot.margin=unit(c(0.2,0.1,0.1,0.1), "cm")
        ,axis.text=element_text(size=6)
        ,legend.key.size = unit(0.5,"line")  
        ,legend.background = element_blank()
        ,legend.text=element_text(size=8)
        ,strip.text=element_text(size=8)
        ,strip.background =element_blank()
  )
dev.off()







#***************************effect of parameter values on the difference between dissemination & transmission************
# 
hcRes <- readRDS(here(".//stochasticModelsSimplified//sensitivityHS240823.rds"))   # higher range of escape rate
names(hcRes)[1:6] <- c("sevenDaysH","tenDaysH","fourteenDaysH","sevenDaysS","tenDaysS","fourteenDaysS")

hcRes <- melt(hcRes,id.vars = c("sevenDaysH","tenDaysH","fourteenDaysH","sevenDaysS","tenDaysS","fourteenDaysS","run"))

hcRes$variable <- factor(hcRes$variable
                         , levels = c("bmClear"
                                      ,"probabilityInfection"
                                      ,"growthRate"
                                      ,"carryCap"
                                      ,"escapeRate"
                         )
                         , ordered = TRUE
                         , labels=c("Bloodmeal clearance rate"
                                    ,"Probability of infection"
                                    ,"Growth rate"
                                    ,"Carrying capacity"
                                    ,"Escape rate")
                         
)

# probably better to code this in and do on numerators, but just use proportions for now
hcRes$diffHS7 <- hcRes$sevenDaysH - hcRes$sevenDaysS
hcRes$diffHS10 <- hcRes$tenDaysH - hcRes$tenDaysS
hcRes$diffHS14 <- hcRes$fourteenDaysH - hcRes$fourteenDaysS

tiff(here(".//stochasticModelsSimplified//fig_HSevenDaysSensitivity2.tiff"), height =4, width = 5, units = 'in', compression="lzw", res=400)
ggplot(hcRes) +
  geom_point(aes(x=value,y=sevenDaysH),col="grey",alpha=0.5,size=0.5) +
  facet_wrap(~variable, scales="free") +
  xlab("Parameter value") +
  ylab("Proportion of simulations with hemocoel infection") +
  theme_bw() +
  theme(panel.border = element_blank()                   
        ,axis.line = element_line(color = 'black')
        ,text=element_text(size=10)
        ,plot.margin=unit(c(0.2,0.1,0.1,0.1), "cm")
        ,axis.text=element_text(size=6)
        ,legend.key.size = unit(0.5,"line")  
        ,legend.background = element_blank()
        ,legend.text=element_text(size=8)
        ,strip.text=element_text(size=8)
        ,strip.background =element_blank()
  )
dev.off()


tiff(here(".//stochasticModelsSimplified//fig_StenDaysSensitivity2.tiff"), height =4, width = 5, units = 'in', compression="lzw", res=400)
ggplot(hcRes) +
  geom_point(aes(x=value,y=tenDaysS),col="grey",alpha=0.5,size=0.5) +
  facet_wrap(~variable, scales="free") +
  xlab("Parameter value") +
  ylab("Proportion of simulations with salivary gland infection") +
  theme_bw() +
  theme(panel.border = element_blank()                   
        ,axis.line = element_line(color = 'black')
        ,text=element_text(size=10)
        ,plot.margin=unit(c(0.2,0.1,0.1,0.1), "cm")
        ,axis.text=element_text(size=6)
        ,legend.key.size = unit(0.5,"line")  
        ,legend.background = element_blank()
        ,legend.text=element_text(size=8)
        ,strip.text=element_text(size=8)
        ,strip.background =element_blank()
  )
dev.off()




tiff(here(".//stochasticModelsSimplified//fig_SsevenDaysSensitivityDifference.tiff"), height =4, width = 5, units = 'in', compression="lzw", res=400)
ggplot(hcRes) +
  geom_point(aes(x=value,y=diffHS7),col="grey",alpha=0.5,size=0.5) +
  facet_wrap(~variable, scales="free") +
  xlab("Parameter value") +
  ylab("Proportion disseminated - proportion transmitting") +
  theme_bw() +
  theme(panel.border = element_blank()                   
        ,axis.line = element_line(color = 'black')
        ,text=element_text(size=10)
        ,plot.margin=unit(c(0.2,0.1,0.1,0.1), "cm")
        ,axis.text=element_text(size=6)
        ,legend.key.size = unit(0.5,"line")  
        ,legend.background = element_blank()
        ,legend.text=element_text(size=8)
        ,strip.text=element_text(size=8)
        ,strip.background =element_blank()
  )
dev.off()


