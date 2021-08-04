# ---------------------------------------------------
# dB - Change in BIOMASS
# -----------------------------------------------------

dBdt.s <- function(t, B, params){
  with(as.list(params), { 
    P <- photosyn.shade(T., I., PT.max, PL.max, T.opt, T.max, Ik, Ca = B, K)
    # This is assuming that T. and I. are constant
    #Need to look up how to make them vary. 
    R <- Resp.temp(R.max, T., RT.opt, RT.max)
    dB <- (P - R - M)*B
    list(dB)
  })
}