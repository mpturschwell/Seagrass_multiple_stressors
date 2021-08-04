######################################################################################
# Function for attack rate and mortality scaling with temperate
#####################################################################################
# Lopez-Urrutia 2008 - The metabolic theory of ecology and algal bloom formation
Attack.Mortal.Temp.Scaling <- function(T.){
  # T. is the temperature in celsius
  b0 <- 43069437027 # 1/Attack.Mortal.Temp.Scaling_unscaled(34.9)
  E <- 0.65 # activation energy of heterotroph as per Lopez-Urrutia 2008
  k <- 8.617333262e-5 # Boltzmann constant as per Wikipedia
  b0*exp(-E/(k*(T. + 273.15))) 
}

