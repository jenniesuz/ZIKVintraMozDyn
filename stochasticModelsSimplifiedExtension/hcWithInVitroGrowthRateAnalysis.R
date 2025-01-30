
library(here)
source(here(".//stochasticModelsSimplifiedExtension//InfDissTransRepeatModelFunc.R"))
library(grid)
library(gridExtra)
library(ggplot2)
library(parallel)

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

#**********************************************************************************


#**********************************************************************************

propsSen1$run <- "1: Baseline"
propsSen2$run <- "2: H growth 1/10 of M"
propsSen3$run <- "3: H growth 1/10 of M plus H growth var"
