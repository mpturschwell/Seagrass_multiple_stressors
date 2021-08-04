#######################################################################################
# PHOTOSYNTHESIS MODEL - Combine light and temp into one model 
#######################################################################################
photosyn.shade <- function(T., I., PT.max, PL.max, T.opt, T.max, Ik, Ca, B.max){
  PL.max <- Photo.Temp(T., PT.max, T.opt, T.max)
  Photo.Light.Shade(I., PL.max, Ik, Ca, B.max)
}