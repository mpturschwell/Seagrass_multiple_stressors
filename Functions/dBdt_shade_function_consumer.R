# ---------------------------------------------------
# dB - Change in BIOMASS
# -----------------------------------------------------

dBdt.s.consumer <- function(t, B, params){
  with(as.list(params), { 
    P <- photosyn.shade(T., I., PT.max, PL.max, T.opt, T.max, Ik, Ca = B, B.max)
    R <- Resp.temp(R.max, T., RT.opt, RT.max)
    dB <- (P - R - M - (a*X))*B
    list(dB)
  })
}