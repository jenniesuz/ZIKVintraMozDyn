library(adaptivetau)

infDissTransModel <- function(startingVirus
                           ,bloodmealClearance = 1/3
                           ,propSuccessInf = 10^-3
                           ,growthRate = 0.1
                           ,carryCap = 10^8
                           ,escapeRate = 0.05
){
  
  
  params <- list(bloodmealClearance = bloodmealClearance
                 ,propSuccessInf = propSuccessInf
                 ,growthRate = growthRate 
                 ,carryCap =  carryCap
                 ,escapeRate = escapeRate    
  )
  
  transitions <- list(c(Gv = -1)
                      ,c(Gv = -1, Mv = +1)
                      ,c(Mv = +1)
                      ,c(Mv= -1, Hv = +1)
                      ,c(Hv = +1)
                      ,c(Hv= -1, Sv = +1)
  )
  
  
  lvrates <- function(y,params,t){
    return( c(params$bloodmealClearance*y["Gv"]               # virus is cleared as bloodmeal digested
              ,y["Gv"]*propSuccessInf
              ,y["Mv"]*params$growthRate*(carryCap - y["Mv"])/params$carryCap
              ,params$escapeRate*y["Mv"]
              ,y["Hv"]*params$growthRate*(carryCap - y["Hv"])/params$carryCap
              ,params$escapeRate*y["Hv"]

    )
    
    )
  }
  
  out <- ssa.adaptivetau(c(Gv = round(startingVirus*0.003,0), Mv = 0, Hv = 0, Sv = 0),
                       transitions, lvrates, params, tf=360
                       , tl.params=list(epsilon=0.005)) # time accuracy trade-off - larger values of epsilon less accurate - may overshoot therefore get negative values but quicker
  return(data.frame(out))
}



