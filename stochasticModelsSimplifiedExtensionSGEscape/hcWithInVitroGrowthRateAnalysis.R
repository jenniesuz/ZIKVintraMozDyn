
# in the laboratory results all chimeras produced a disseminated infection across
# all time points but transmission looked different 

library(here)
source(here(".//stochasticModelsSimplifiedExtensionSGEscape//InfDissTransRepeatModelFunc.R"))
library(grid)
library(gridExtra)
library(ggplot2)
library(parallel)

# tf = 360, epilson = 0.005, starting virus = 10^6

#*********************
m <- 0.001
v <-  0.0000001
sc <- v/m
sh <- m^2/v
rate <- 1/sc
x.max = qgamma(0.999, shape=sh, scale=sc)
x = seq(from=0, to=x.max, by=x.max/1000)
dens = dgamma(x, shape=sh, scale=sc)
plot(x, dens, type='l')


#***************Senegal median growth rate in vitro***********************
#*****    ***
modelOutput <- repeatModel(virus_params(bloodmealClearance = 1/72
                                        ,propSuccessInf = 10^-4
                                        ,growthRateH = 0.04
                                        ,growthRateM = 0.04
                                        ,carryCap = 10^20
                                        ,escapeRateM = 0.00005 
                                        ,escapeRateH = 0.00005/10
                                        ,escapeVarS = 0.000000002
                                        ,escapeVarSTrue = F) # 0.001, 0.00005
                           ,startingVirus=10^6)

modelOutput <- do.call(rbind,modelOutput)

propsSen4 <- dissSummaryFunc(modelOutput)

# variation in the salivary gland infection barrier between simulations & lower overall average could explain
#******************************************************************************
modelOutput <- repeatModel(virus_params(bloodmealClearance = 1/72
                                        ,propSuccessInf = 10^-4
                                        ,growthRateH = 0.04
                                        ,growthRateM = 0.04
                                        ,carryCap = 10^20
                                        ,escapeRateM = 0.00005
                                        ,escapeRateH = 0.0005/10
                                        ,escapeVarS = 0.0000002
                                        ,escapeVarSTrue = T) # 0.001, 0.00005
                           ,startingVirus=10^6)

modelOutput <- do.call(rbind,modelOutput)

propsSen5 <- dissSummaryFunc(modelOutput)

propsSen4$run <- "4: S escape 1/10 of M:H escape"
propsSen5$run <- "5: S escape 1/10 of M:H escape plus var"


#****************************Plot all**************************************************
#*


cols <- c("#8c510a"
          ,"#dfc27d"
          ,"#c7eae5"
          ,"#35978f"
          ,"#01665e"
)


propsSen <- rbind.data.frame(propsSen1,propsSen2,propsSen3,propsSen4,propsSen5)


tiff(here::here(".//stochasticModelsSimplifiedExtensionSGEscape//fig_InfUpdated.tiff")
     , height =6, width = 4, units = 'in', compression="lzw", res=400)


ggplot(propsSen) +
  geom_line(aes(x=time,y=proportionInfected,col=run,group=run),linetype="longdash") +
  xlim(1,14) +
  scale_color_manual(name="run",values=cols) +
  xlab("Time (days)") +
  ylab("Proportion \n infected") +
  theme_set(theme_bw())  +
  theme(panel.border = element_blank()
        ,axis.line = element_line(color = 'black')
        ,text=element_text(size=12)
        ,plot.margin=unit(c(0.2,0.3,0.1,0.1), "cm")
        ,axis.text=element_text(size=10)
        ,legend.key.size = unit(0.5,"line")
        ,legend.background = element_blank()
        ,legend.text=element_text(size=10)
        ,legend.position =c(0.6,0.5)
        ,legend.title = element_blank()
  )



dev.off()

tiff(here::here(".//stochasticModelsSimplifiedExtensionSGEscape//fig_dissUpdated.tiff")
     , height =6, width = 4, units = 'in', compression="lzw", res=400)


ggplot(propsSen) +
  geom_line(aes(x=time,y=proportionDisseminated,col=run,group=run),linetype="longdash") +
  xlim(1,14) +
  scale_color_manual(name="run",values=cols) +
  xlab("Time (days)") +
  ylab("Proportion \n disseminated") +
  theme_set(theme_bw())  +
  theme(panel.border = element_blank()
        ,axis.line = element_line(color = 'black')
        ,text=element_text(size=12)
        ,plot.margin=unit(c(0.2,0.3,0.1,0.1), "cm")
        ,axis.text=element_text(size=10)
        ,legend.key.size = unit(0.5,"line")
        ,legend.background = element_blank()
        ,legend.text=element_text(size=10)
        ,legend.position ="none"
        ,legend.title = element_blank()
  )


dev.off()




tiff(here::here(".//stochasticModelsSimplifiedExtensionSGEscape//fig_transUpdated.tiff")
     , height =6, width = 4, units = 'in', compression="lzw", res=400)


ggplot(propsSen) +
  geom_line(aes(x=time,y=proportionCapableofTransmission,col=run,group=run),linetype="longdash") +
  xlim(1,14) +
  scale_color_manual(name="run",values=cols) +
  xlab("Time (days)") +
  ylab("Proportion with virus \n in salivary glands") +
  theme_set(theme_bw())  +
  theme(panel.border = element_blank()
        ,axis.line = element_line(color = 'black')
        ,text=element_text(size=12)
        ,plot.margin=unit(c(0.2,0.3,0.1,0.1), "cm")
        ,axis.text=element_text(size=10)
        ,legend.key.size = unit(0.5,"line")
        ,legend.background = element_blank()
        ,legend.text=element_text(size=10)
        ,legend.position ="none"
        ,legend.title = element_blank()
  )

#saveRDS(dissPlotHemo,"dissPlotHemo.rds")
#saveRDS(transPlotHemo,"transPlotHemo.rds")



dev.off()




