# ---------------------------------------------------
# dB - Change in BIOMASS
# -----------------------------------------------------

dBdt.s.consumer.foodweb <- function(t, B, X, params){
  # Base model
  with(as.list(params), { 
    P <- photosyn.shade(T., I., PT.max, PL.max, T.opt, T.max, Ik, Ca = B, B.max)
    R <- Resp.temp(R.max, T., RT.opt, RT.max)
    S <- Attack.Mortal.Temp.Scaling(T.)
    dB <- (P - R - M - (a*S*X))*B
    dX <- (a*S*c*X*B)-(v*S*X)
    c(dB,dX)
  })
}

