library(here)
library(ggplot2)

dat <- read.csv(here(".//data//growthCurveForR.csv"))

head(dat)

datLong <- reshape2::melt(dat,id.var=c("chimera","replicate"))
names(datLong) <- c("Chimera","Replicate","Time","PFU")
temp <- strsplit(as.character(datLong$Time),"X")
temp <- do.call(rbind,temp)
temp <- temp[,2]
datLong$Time <- as.numeric(temp)
head(datLong)



cols <- c("#7fc97f"
  ,"#beaed4"
  ,"#fdc086"
  ,"#ffff99"
  ,"#386cb0"
  ,"#f0027f"
  ,"#bf5b17"
  ,"#666666"
)

dataPlot <- ggplot(datLong) +
  geom_line(aes(x=Time 
                ,y=log10(PFU)
                ,col=Chimera
                ,group=interaction(Chimera,Replicate))) +
  xlim(24,50) +
  scale_color_manual("Chimera", values = cols) +
  theme_grey()

dataPlot