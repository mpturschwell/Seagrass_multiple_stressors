#----------------------------------------------
# constant to convert umol Oxygen to mg carbon 
#----------------------------------------------

# umol to mg = *0.000001
# multiply by molecular weight of oxygen (O2 = 32g)
# assume photosynthetic quotient of 1: 1 mole O2 = 1 mole C
# = weight of 1 mole CO2(44g) / 1 mole O2(32g))
umol_to_mg <- function(umol) (umol*0.000001 * 32)*1.375  


#----------------------------------------------
# constant to convert R.max and   
#----------------------------------------------

# Conversion factor for:   mg C h-1  ->  g dry wt day-1
# MISCHA add details here
mg.C.per.h_to_g.dry.wt.per.d <- function(x) x * 24 * (1 / 0.33) / 1000 