library(adaptivetau)

infDissTransModel <- function(startingVirus
                           ,bloodmealClearance = 1/3
                           ,propSuccessInf = 10^-3
                           ,growthRateM = 0.1
                           ,growthRateH = 0.05
                           ,carryCap = 10^8
                           ,escapeRateM = 0.05
                           ,escapeRateH = 0.05
){
  
  
  params <- list(bloodmealClearance = bloodmealClearance
                 ,propSuccessInf = propSuccessInf
                 ,growthRateM = growthRateM
                 ,growthRateH = growthRateH
                 ,carryCap =  carryCap
                 ,escapeRateM = escapeRateM
                 ,escapeRateH = escapeRateH
  )
  
  transitions <- list(c(Gv = -1)
                      ,c(Gv = -1, Mv = +1)
                      ,c(Mv = +1)
                      ,c(Hv = +1)
                      ,c(Hv = +1)
                      ,c(Sv = +1)
  )
  
  
  lvrates <- function(y,params,t){
    return( c(params$bloodmealClearance*y["Gv"]              
              ,y["Gv"]*propSuccessInf
              ,y["Mv"]*params$growthRateM*(carryCap - y["Mv"])/params$carryCap
              ,params$escapeRateM*y["Mv"]
              ,y["Hv"]*params$growthRateH*(carryCap - y["Hv"])/params$carryCap
              ,params$escapeRateH*y["Hv"]

    )
    
    )
  }
  
  out <- ssa.adaptivetau(c(Gv = round(startingVirus*0.003,0), Mv = 0, Hv = 0, Sv = 0),
                       transitions, lvrates, params, tf=360
                       , tl.params=list(epsilon=0.005)) 
  return(data.frame(out))
}



