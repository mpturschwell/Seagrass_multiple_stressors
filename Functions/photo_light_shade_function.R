#######################################################################################
# GROWTH - PHOTOSYNTHESIS - LIGHT function with self shading
# Reference: Burd & Dunton 2001 - MEPS
# Species: Halodule wrightii
#######################################################################################
#' Jassby-Platt parameterization (Jassby & Platt 1976): gross production (P) mmol O2 g dry wt-1 h-1 at irradiance (I) mmol m-2 s-1
#' @param I. 
#' 
#' 
Photo.Light.Shade <- function(I., PL.max, Ik, Ca, B.max){ 
  PL.max*tanh(I./Ik)*(1-Ca/B.max) 
}


