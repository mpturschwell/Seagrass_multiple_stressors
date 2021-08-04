# ---------------------------------------------------
# dB - Change in BIOMASS
# -----------------------------------------------------

dBdt.s.consumer.foodweb <- function(t, B, X, params){
  with(as.list(params), { 
    P <- photosyn.shade(T., I., PT.max, PL.max, T.opt, T.max, Ik, Ca = B, B.max)
    R <- Resp.temp(R.max, T., RT.opt, RT.max)
    S <- Attack.Mortal.Temp.Scaling(T.)
    #dB <- ((P - R - M - (a*S*X))*B)*24*(1/0.33)/1000 # convert P & R (mg/hr) to g (/1000) of seagrass tissue (x3) per day (*24)
    dB <- (P - R - M - (a*S*X))*B
    #dB <- ((P - R)*24*(1/0.33)/1000 - M - (a*S*X))*B 
    dX <- (a*S*c*X*B)-(v*S*X)
    c(dB,dX)
  })
}




#dBdt.s.consumer.foodweb.X <- function(t, B, X, params){
#  with(as.list(params), { 
#    dX <- (a*c*X*B)-(v*X)
#    list(dX)
#  })
# }
