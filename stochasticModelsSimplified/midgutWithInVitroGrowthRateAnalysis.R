
library(here)
source(here(".//stochasticModelsSimplified//InfDissTransRepeatModelFunc.R"))
library(grid)
library(gridExtra)
library(ggplot2)
library(parallel)

# tf = 168, epilson = 0.005, starting virus = 10^6, 100 mosquitoes per sim

source(here(".//data//r_plotDoseResponse.R"))

infPrev <- infPrev[infPrev$Chimera %in% c("Senegal","Thai"),]

obsInfPlot <- ggplot(infPrev) +
  geom_point(aes(x=log10(Dose),y=mean,shape=Chimera),col="grey") +
  geom_errorbar(aes(x=log10(Dose),ymin=lower,ymax=upper),col="grey") +
  xlim(2,7) +
  ylab("Proportion of mosquitoes infected") +
  theme_set(theme_bw())  +
  theme(panel.border = element_blank()
        ,axis.line = element_line(color = 'black')
        ,text=element_text(size=10)
        ,plot.margin=unit(c(0.2,0.3,0.1,0.1), "cm")
        ,axis.text=element_text(size=10)
        ,legend.key.size = unit(0.8,"line")
        ,legend.background = element_blank()
        ,legend.text=element_text(size=10)
        #    ,legend.position =c(9,0.5)
        #  ,legend.title = element_blank()
  )


obsInfPlot



#**************************************
#*****************Beta 10^-3.5**************************************

virConcs <- c(10^2,10^3,10^3.4,10^3.8,10^4,10^4.4,10^4.8,10^5,10^5.4,10^5.8,10^6,10^6.5,10^7)

cl <- makeCluster(detectCores()-1)               # get cores from your computer - 1
clusterEvalQ(cl, {library(adaptivetau)})
environment(repeatModel) <- .GlobalEnv        # need to figure out whether this is needed
clusterExport(cl, varlist=c("infDissTransModel"
                            ,"virus_params"
                            ,"repeatModel"
                            ,"virConcs"
                  
),
envir=environment())               # each cluster doesn't have access to any of the code we have run above
# so we need to send it a list of objects 

start <- Sys.time()
simsMidgut <- parLapply(cl,1:length(virConcs),function(y){    # this is lapply but in parallel - so it runs one simulation on each cluster
  
  vir <- virConcs[y]                 # take the yth row
  
  modelOutput <- repeatModel(virus_params(bloodmealClearance = 1/72
                                          ,propSuccessInf = 10^-3.5
                                          ,growthRate = 0.04
                                          ,carryCap = 10^20
                                          ,escapeRate = 0.00005)
                             ,startingVirus=vir)
  
  modelOutput <- do.call(rbind,modelOutput)
  # midgut proprtion infection
  propInf <- sum(modelOutput$inf[!duplicated(modelOutput$run)])/length(modelOutput$inf[!duplicated(modelOutput$run)])
  inf <- sum(modelOutput$inf[!duplicated(modelOutput$run)])
  den <- length(modelOutput$inf[!duplicated(modelOutput$run)])
  return(c(vir,propInf,inf,den) )                            # return
})
stopCluster(cl)   # stop the clusters
end <- Sys.time()
end - start

simsMidgut2 <- do.call(rbind.data.frame,simsMidgut)
names(simsMidgut2) <- c("conc","propInf","inf","den")
dr10tom3 <- simsMidgut2
#**********************************************************************


#*****************Beta 10^-4.5**************************************
virConcs <- c(10^2,10^3,10^3.4,10^3.8,10^4,10^4.4,10^4.8,10^5,10^5.4,10^5.8,10^6,10^6.5,10^7)

cl <- makeCluster(detectCores()-1)               # get cores from your computer - 1
clusterEvalQ(cl, {library(adaptivetau)})
environment(repeatModel) <- .GlobalEnv        # need to figure out whether this is needed
clusterExport(cl, varlist=c("infDissTransModel"
                            ,"virus_params"
                            ,"repeatModel"
                            ,"virConcs"
),
envir=environment())               # each cluster doesn't have access to any of the code we have run above
# so we need to send it a list of objects 

start <- Sys.time()
simsMidgut <- parLapply(cl,1:length(virConcs),function(y){    # this is lapply but in parallel - so it runs one simulation on each cluster
  
  vir <- virConcs[y]                 # take the yth row
  
  modelOutput <- repeatModel(virus_params(bloodmealClearance = 1/72
                                          ,propSuccessInf = 10^-4.25
                                          ,growthRate = 0.04
                                          ,carryCap = 10^20
                                          ,escapeRate = 0.00005)
                             ,startingVirus=vir)
  
  modelOutput <- do.call(rbind,modelOutput)
  # midgut proprtion infection
  propInf <- sum(modelOutput$inf[!duplicated(modelOutput$run)])/length(modelOutput$inf[!duplicated(modelOutput$run)])
  inf <- sum(modelOutput$inf[!duplicated(modelOutput$run)])
  den <- length(modelOutput$inf[!duplicated(modelOutput$run)])
  return(c(vir,propInf,inf,den) )                            # return
})
stopCluster(cl)   # stop the clusters
end <- Sys.time()
end - start

simsMidgut2 <- do.call(rbind.data.frame,simsMidgut)
names(simsMidgut2) <- c("conc","propInf","inf","den")
dr10tom4 <- simsMidgut2
#**********************************************************************



#*****************Beta 10^-5**************************************
virConcs <- c(10^2,10^3,10^3.4,10^3.8,10^4,10^4.4,10^4.8,10^5,10^5.4,10^5.8,10^6,10^6.5,10^7)

cl <- makeCluster(detectCores()-1)               # get cores from your computer - 1
clusterEvalQ(cl, {library(adaptivetau)})
environment(repeatModel) <- .GlobalEnv        # need to figure out whether this is needed
clusterExport(cl, varlist=c("infDissTransModel"
                            ,"virus_params"
                            ,"repeatModel"
                            ,"virConcs"
),
envir=environment())               # each cluster doesn't have access to any of the code we have run above
# so we need to send it a list of objects 

start <- Sys.time()
simsMidgut <- parLapply(cl,1:length(virConcs),function(y){    # this is lapply but in parallel - so it runs one simulation on each cluster
  
  vir <- virConcs[y]                 # take the yth row
  
  modelOutput <- repeatModel(virus_params(bloodmealClearance = 1/72
                                          ,propSuccessInf = 10^-5
                                          ,growthRate = 0.04
                                          ,carryCap = 10^20
                                          ,escapeRate = 0.00005)
                             ,startingVirus=vir)
  
  modelOutput <- do.call(rbind,modelOutput)
  # midgut proprtion infection
  propInf <- sum(modelOutput$inf[!duplicated(modelOutput$run)])/length(modelOutput$inf[!duplicated(modelOutput$run)])
  inf <- sum(modelOutput$inf[!duplicated(modelOutput$run)])
  den <- length(modelOutput$inf[!duplicated(modelOutput$run)])
  return(c(vir,propInf,inf,den) )                            # return
})
stopCluster(cl)   # stop the clusters
end <- Sys.time()
end - start

simsMidgut2 <- do.call(rbind.data.frame,simsMidgut)
names(simsMidgut2) <- c("conc","propInf","inf","den")
dr10tom5 <- simsMidgut2
#**********************************************************************

#EFA5A5

dr10tom3$Beta <- "-3.5"
dr10tom4$Beta <- "-4.25"
dr10tom5$Beta <- "-5"
mod <- rbind.data.frame(dr10tom3,dr10tom4,dr10tom5)

tiff(here::here(".//stochasticModelsSimplified//fig_midgutDoseResponse.tiff")
     , height =4, width = 5, units = 'in', compression="lzw", res=400)

obsInfPlot +
  geom_line(data=mod,aes(x=log10(conc),y=propInf,col=Beta),linetype="longdash",size=1.1,alpha=1) +
  scale_color_manual(values=c("#8c510a","#d8b365","#01665e")) +
  ylab("Proportion of simulations with midgut infection") +
  xlab("Virus concentration (log 10)") +
  labs(col="Beta (log10)") +
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

# Sen 10^-4
# Tha 10^-4.5

# growth rate doesnt matter
# some evidence of small effect of nsp and or sp - stats anal of data