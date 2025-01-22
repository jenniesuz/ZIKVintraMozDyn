library(here)
library(ggplot2 )
library(reshape2)



#*******************************Hemocoel & Salivary glad with In Vitro Growth Rate***********************
#*
hcResSenegal <- readRDS(here(".//stochasticModelsSimplified//HSSenegalInvitroGrowthRatees0.0005.rds"))
hcResThai <- readRDS(here(".//stochasticModelsSimplified//HSThaiInvitroGrowthRatees0.0005.rds"))

names(hcResSenegal)<- c("sevenDaysH","tenDaysH","fourteenDaysH","sevenDaysS","tenDaysS","fourteenDaysS","GrowthRate")
names(hcResThai)<- c("sevenDaysH","tenDaysH","fourteenDaysH","sevenDaysS","tenDaysS","fourteenDaysS","GrowthRate")

hcResSenegal <- melt(hcResSenegal,id.vars = c("sevenDaysH","tenDaysH","fourteenDaysH","sevenDaysS","tenDaysS","fourteenDaysS"))
hcResThai <- melt(hcResThai,id.vars = c("sevenDaysH","tenDaysH","fourteenDaysH","sevenDaysS","tenDaysS","fourteenDaysS"))

hcResSenegal$Chimera <- "Senegal"
hcResThai$Chimera <- "Thai"
hcRes <- rbind.data.frame(hcResSenegal,hcResThai)

#***************************Plot effect of growth rate on diss/ trans***********************
ggplot(hcRes) +
  geom_point(aes(x=value,y=sevenDaysH),col="grey",alpha=0.5,size=0.5) +
  geom_point(aes(x=value,y=sevenDaysS),col="red",alpha=0.5,size=0.5) +
  facet_wrap(~variable, scales="free") +
  xlab("Parameter value") +
  ylab("Proportion of simulations with dissemination/ transmission") +
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
  ) +
  facet_wrap(~Chimera)



#tiff(here(".//stochasticModelsSimplified//fig_HSevenDaysSensitivity.tiff"), height =4, width = 5, units = 'in', compression="lzw", res=400)
ggplot(hcRes) +
  geom_histogram(aes(x=sevenDaysH,fill=Chimera),alpha=0.5)

ggplot(hcRes) +
  geom_histogram(aes(x=sevenDaysS,fill=Chimera),alpha=0.5)



#***********************************
hcResSenegal <- readRDS(here(".//stochasticModelsSimplified//HSSenegalInvitroGrowthRatees0.0005.rds"))
hcResThai <- readRDS(here(".//stochasticModelsSimplified//HSThaiInvitroGrowthRatees0.0005.rds"))

names(hcResSenegal)<- c("sevenDaysH","tenDaysH","fourteenDaysH","sevenDaysS","tenDaysS","fourteenDaysS","GrowthRate")
names(hcResThai)<- c("sevenDaysH","tenDaysH","fourteenDaysH","sevenDaysS","tenDaysS","fourteenDaysS","GrowthRate")

hcResSenegal <- melt(hcResSenegal,id.vars = c("GrowthRate" ))
hcResThai <- melt(hcResThai,id.vars = c("GrowthRate" ))

hcResSenegal$Chimera <- "Senegal"
hcResThai$Chimera <- "Thai"
hcRes <- rbind.data.frame(hcResSenegal,hcResThai)

tissue <- strsplit(as.character(hcRes$variable),split="ys")
tissue <- sapply(tissue, "[[", 2)

hcRes$tissue <- tissue
hcRes$day <- NA
hcRes$day[hcRes$variable %in% "sevenDaysH"] <- 7
hcRes$day[hcRes$variable %in% "sevenDaysS"] <- 7
hcRes$day[hcRes$variable %in% "tenDaysH"] <- 10
hcRes$day[hcRes$variable %in% "tenDaysS"] <- 10
hcRes$day[hcRes$variable %in% "fourteenDaysH"] <- 14
hcRes$day[hcRes$variable %in% "fourteenDaysS"] <- 14

ggplot(hcRes) +
  geom_point(aes(x=day,y=value,col=tissue)) +
  facet_wrap(~Chimera)

# need a way to summarise all the sims
# take numerator from the output rather than the proportion?
library(plyr)
resProp <- ddply(hcRes,.(Chimera,tissue,day),summarise,mp=mean(value))
ggplot(resProp) +
  geom_point(aes(x=day,y=mp,col=Chimera)) +
  xlim(0,14) +
  facet_wrap(~tissue)

# for now do fewer sampling from the posterior
# then simulate and save the results from more days?