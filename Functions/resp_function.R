####################################################################################################
# RESPIRATION - TEMPERATURE function 
# can be represented by the same P.T function according to Collier et al. 2017 - C. serrulata
####################################################################################################
Resp.temp <- function(R.max, T., RT.opt, RT.max){ 
  R.max*((RT.max-T.)/(RT.max-RT.opt)) * ((T./RT.opt) ^ (RT.opt/(RT.max-RT.opt))) 
}
