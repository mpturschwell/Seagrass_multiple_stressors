# ---------------------------------
# Biomass integration 
# ---------------------------------
# looks at chnage in biomass through time based on netgrowth and mortality 
# tmax - length of vector, dt - timestep, B1 = intial biomass
# calling dBdt.s function with calculates chnage in biomass based on net growth - mortality 
biomass.int.s.consumer <- function(tmax, dt, B_init, params){
  
  B <- rep(NA, tmax) # create empyt vector to store output
  B[1] <- B_init         # set inital biomass in time 1
  
  for (time_step in 1:(tmax-1)){
    B[time_step+1] <- B[time_step] + 
      dBdt.s.consumer(1, B[time_step], params)[[1]]*dt
  }
  return(B)
}

# create wrapper to list parameters before calling pmap 
biomass.int.wrapper.cons <- function(tmax, dt, B_init, T., I., PT.max, PL.max, T.opt, T.max, Ik, B.max, 
                                R.max, RT.opt, RT.max, M, a, X){
  params <- list(
      B.max = B.max, 
      T. = T.,     
      I. = I.,     
      PT.max = PT.max,
      PL.max = PL.max,
      T.opt = T.opt,
      T.max = T.max,
      Ik = Ik,
      R.max = R.max,
      RT.opt = RT.opt,
      RT.max = RT.max,  
      M = M,
      a = a,
      X = X
    )
  
  biomass.int.s.consumer(tmax, dt, B_init, params)
  
}
