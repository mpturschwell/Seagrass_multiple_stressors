# ---------------------------------------------------
# dB - Change in BIOMASS
# -----------------------------------------------------

dBdt.s.consumer <- function(t, B, params){
  with(as.list(params), { 
    P <- photosyn.shade(T., I., PT.max, PL.max, T.opt, T.max, Ik, Ca = B, K)
    R <- Resp.temp(R.max, T., RT.opt, RT.max)
    dB <- ((P - R - M - (a*X))*B)*24*(1/0.33)/1000 # convett to g (/1000) seagrass tissue (x3) per day (*24) 
    list(dB)
  })
}