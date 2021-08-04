######################################################################################
# GROWTH - TEMPERATURE function: Adams et al. 2017 - Scientific Reports: sp. C. serrulata
#####################################################################################
# Yan and Hunt model (Yan & Hunt 1999): gross photosynthesis (P) in mg C g-1 DW h-1 at time T
Photo.Temp <- function(T., PT.max, T.opt, T.max){ 
  PT.max*((T.max-T.)/(T.max-T.opt)) * ((T./T.opt) ^ (T.opt/(T.max-T.opt))) 
}
