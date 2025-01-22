
# in the laboratory results all chimeras produced a disseminated infection across
# all time points but transmission looked different 

library(here)
source(here(".//stochasticModelsSimplifiedExtension//InfDissTransRepeatModelFunc.R"))
library(grid)
library(gridExtra)
library(ggplot2)
library(parallel)


#***************Just take Senegal growth rate in vitro***********************
#*
#*****no between-mosquito variation in hemocoel growth rate and same growth rate across tissues***
modelOutput <- repeatModel(virus_params(bloodmealClearance = 1/72
                                        ,propSuccessInf = 10^-4
                                        ,growthRateM = 0.04
                                        ,growthRateH = 0.04
                                        ,carryCap = 10^20
                                        ,escapeRateM = 0.00005 #0.0005
                                        ,escapeRateH = 0.00005 #0.0005
                                        ,growthVarH = 0
                                        ,growthVarHTrue = F) # 0.001, 0.00005
                           ,startingVirus=10^6)

modelOutput <- do.call(rbind,modelOutput)


propsSen1 <- dissSummaryFunc(modelOutput)
#******************************************************************************


#********no between-mosquito variation in hemocoel growth rate but lower growth rate in hemocoeol**********

#************************************
# 1/10th
modelOutput <- repeatModel(virus_params(bloodmealClearance = 1/72
                                        ,propSuccessInf = 10^-4
                                        ,growthRateM = 0.04
                                        ,growthRateH = 0.04/10
                                        ,carryCap = 10^20
                                        ,escapeRateM =  0.00005
                                        ,escapeRateH =  0.00005
                                        ,growthVarH = 0
                                        ,growthVarHTrue = F) # 0.001, 0.00005
                           ,startingVirus=10^6)

modelOutput <- do.call(rbind,modelOutput)

propsSen2 <- dissSummaryFunc(modelOutput)

# tends to shift both curves to the right - because hemocoel growth rate also affects dissemination
#**********************************************************************************




#********no between-mosquito variation, same growth rate but different escape rates**********


#1/10th
modelOutput <- repeatModel(virus_params(bloodmealClearance = 1/72
                                        ,propSuccessInf = 10^-4
                                        ,growthRateM = 0.04
                                        ,growthRateH = 0.04/10
                                        ,carryCap = 10^20
                                        ,escapeRateM =  0.00005
                                        ,escapeRateH =  0.00005
                                        ,growthVarH = 0.00001
                                        ,growthVarHTrue = T) # 0.001, 0.00005
                           ,startingVirus=10^6)

modelOutput <- do.call(rbind,modelOutput)

propsSen3 <- dissSummaryFunc(modelOutput)



# lower escape rate to SG just shifts SG curve further to the right but doesn't affect the slope of the curve
#**********************************************************************************


  
  # This can flatten the curve, but as above, because it also affects dissemination, then it affects both H and SG curves
  # and it requires a lot of variation to create this difference?
  
  # makes the lab findings of the difference between dissemination and SG infection quite astounding?
  
  # did we need a model to be able to come to this conclusion though?! ...
  # dose response and dissemination conform to the very simple model process of stochastcity without additional
  # complexities, but the process after that, for virus to get to the saliva is more complicated
  
  # essentially need in the model some mechanism that prevents some mosquitoes that have a disseminated infection
  # from ever getting a salivary gland infection.
  #**********************************************************************************

  

propsSen1$run <- "1: Baseline"
propsSen2$run <- "2: H growth 1/10 of M"
propsSen3$run <- "3: H growth 1/10 of M plus H growth var"
