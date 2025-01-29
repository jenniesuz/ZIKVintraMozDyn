library("here")

source(here(".//inVitro//dataFirstSet.R"))

dat <- dat[dat$Chimera %in% c("Senegal","Thai"),]

pGrowth <- ggplot(dat) +
  geom_point(aes(x=time
                 ,y=log10(Titer)
                 ,shape=Chimera
  ),col="grey") +
  xlim(0,100) +
  ylim(0,10) +
  theme_bw()


growthFunc <- function(time,k,s,r){
  virus <- s*k*exp(r*time)/((k-s)+s*exp(r*time))
  return(virus)
}

growthFuncV <- Vectorize(growthFunc)


virusGrowth <- growthFunc(time=c(1:100),r=0.035,k=19.57197,s=2.302585)
vg1 <- cbind.data.frame(virusGrowth,time=1:100)

virusGrowth <- growthFunc(time=c(1:100),r=0.05,k=19.57197,s=2.993361)
vg2 <- cbind.data.frame(virusGrowth,time=1:100)

simDat <- cbind.data.frame(time=rep(1:100,2)
                           ,titer=c(vg1$virusGrowth,vg2$virusGrowth)
                           ,parms=c(rep("r = 0.035, s = 1",100)
                                    ,rep("r = 0.05, s = 1.3",100)))


log10(exp(20))
log10(exp(2))
log10(exp(3))

log(10^8.5)
log(10^1)
log(10^1.3)


cols<- c("#d8b365","#8c510a")


tiff(here::here(".//inVitro//JensUpdateLogisticFunction//fig_logisticModel.tiff")
     , height =4, width = 5, units = 'in', compression="lzw", res=400)
pGrowth +
  geom_line(data=simDat,aes(x=time,y=log10(exp(titer)),col=parms,group=parms),size=1.1) +
  ylab(expression(paste("Virus titer (log"[10]," PFU/ml)"))) +
  xlab("Time (hours)") +
  scale_color_manual(name="Parameter values",values=cols) +
  theme_set(theme_bw())  +
  theme(panel.border = element_blank()
        ,axis.line = element_line(color = 'black')
        ,text=element_text(size=10)
        ,plot.margin=unit(c(0.2,0.3,0.1,0.1), "cm")
        ,axis.text=element_text(size=10)
        ,legend.key.size = unit(0.8,"line")
        ,legend.background = element_blank()
        ,legend.text=element_text(size=10)
        ,legend.position =c(0.2,0.8)
        #  ,legend.title = element_blank()
  )

dev.off()
