#----------------------------------------------
# constant to convert umol Oxygen to mg carbon 
#----------------------------------------------

# umol to mg = *0.000001
# multiply by molecular weight of oxygen (O2 = 32g)
# assume photosynthetic quotient of 1: 1 mole O2 = 1 mole C

umol_to_mg <- function(umol){
  (umol*0.000001 * 32)*1.375  # = weight of 1 mole CO2(44g) / 1 mole O2(32g))
}
