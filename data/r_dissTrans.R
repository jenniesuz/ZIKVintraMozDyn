library(here)
library(binom)

dat <- read.csv(here(".\\data\\030823_Chimera_1-8 time to dissemination.csv"))
head(dat)

dat2 <- lapply(unique(dat$Virus),function(x){
  temp <- dat[dat$Virus %in% x,]
  seven <- temp[temp$Day %in% 7,]
  totSeven <- colSums(seven[5:8])
  sampSeven <- length(seven$Infection)
  
  ten <- temp[temp$Day %in% 10,]
  totTen<- colSums(ten[5:8])
  sampTen <- length(ten$Infection)
  
  fourt <- temp[temp$Day %in% 14,]
  totFourt<- colSums(fourt[5:8])
  sampFourt <- length(fourt$Infection)
  
  Virus <- rep(x,3)
  Samp <- c(sampSeven,sampTen,sampFourt)
  Infection <- c(totSeven[1],totTen[1],totFourt[1])
  Dissemination <- c(totSeven[2],totTen[2],totFourt[2])
  Transmission <- c(totSeven[3],totTen[3],totFourt[3])
  day <- c(7,10,14)
  return(cbind.data.frame(Virus,Samp,Infection,Dissemination,Transmission,day))
  
})

dat2 <- do.call(rbind,dat2)

ggplot(dat2) +
  geom_point(aes(x=day,y=Transmission/Dissemination)) +
  facet_wrap(~Virus) +
  ylim(0,1)
