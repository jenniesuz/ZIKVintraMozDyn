
library(here)
source(here(".//stochasticModelsSimplifiedExtension//InfDissTransRepeatModelFunc.R"))
library(grid)
library(gridExtra)
library(ggplot2)

virus_params <- function(   bloodmealClearance = 1/3
                            ,propSuccessInf = 10^-3
                            ,growthRate = 0.01
                            ,carryCap = 10^10
                            ,escapeRate = 0.05               
)
  return(as.list(environment()))

#***************test code**************


# test run at one single input virus amount (Gv):
startTime <- Sys.time()
test <- repeatModel(virus_params(bloodmealClearance = 1/24
                                 ,propSuccessInf = 10^-5
                                 ,growthRate = 0.05
                                 ,carryCap = 10^8
                                 ,escapeRate = 0.001)
                    ,startingVirus=10^6.2)
endTime <- Sys.time()
endTime - startTime

test <- do.call(rbind,test)

head(test)

# midgut proprtion infection
sum(test$inf[!duplicated(test$run)])/length(test$inf[!duplicated(test$run)])


#************Work out proportions*************

#Work out proportions:
testProps <- dissSummaryFunc(test)

testPropsLong <- cbind.data.frame("prop"=c(testProps$proportionDisseminated
                                           ,testProps$proportionCapableofTransmission)
                                  ,"condition"=c(rep("Proportion \n disseminated",length(testProps[,1]))
                                                 ,rep("Proportion capable \n of transmitting",length(testProps[,1])))
                                  ,"time"=rep(testProps$time,2)
)

ggplot(testPropsLong) +
  geom_line(aes(x=time,y=prop,col=condition)) +
  xlab("Time") +
  ylab("Proportion of simulations") +
  ylim(0,1)




#plot dynamics over time

ggplot(test) +
  geom_line(aes(x=time/24,y=log10(Gv),group=run))

ggplot(test) +
  geom_line(aes(x=time/24,y=log10(Mv),group=run))

ggplot(test) +
  geom_line(aes(x=time/24,y=log10(Hv),group=run))

ggplot(test) +
  geom_line(aes(x=time/24,y=log10(Sv),group=run))


