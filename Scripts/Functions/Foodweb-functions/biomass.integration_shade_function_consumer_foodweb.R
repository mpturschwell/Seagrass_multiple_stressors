# ---------------------------------
# Biomass integration 
# ---------------------------------
# looks at change in biomass through time based on netgrowth and mortality 
# tmax - length of vector
# dt - timestep
# B1 = initial biomass

# calling dBdt.s function with calculates change in biomass based on net growth - mortality 
biomass.int.s.consumer.foodweb <- function(tmax, dt, B_init, X_init, params){
          B <- rep(NA, tmax) # create empyt vector to store output
          B[1] <- B_init         # set inital biomass in time 1
  
          X <- rep(NA, tmax) # create empyt vector to store output
          X[1] <- X_init 
  
  for (time_step in 1:(tmax-1)){
        dBdX <- dBdt.s.consumer.foodweb(1, B[time_step], X[time_step], params)
          B[time_step+1] <- B[time_step] + dBdX[1]*dt
          X[time_step+1] <- X[time_step] + dBdX[2]*dt
      }
 mylist <- list(B,X)
 return(mylist)
}

# create wrapper to list parameters before calling pmap 
biomass.int.wrapper.cons.foodweb <- function(tmax, dt, B_init, X_init, T., I., PT.max, PL.max, T.opt, T.max, Ik, B.max, 
                                R.max, RT.opt, RT.max, M, a, c, v){ 
  params <- list(
      T. = T.,     
      I. = I.,     
      PT.max = PT.max,
      PL.max = PL.max,
      T.opt = T.opt,
      T.max = T.max,
      Ik = Ik,
      B.max = B.max, 
      R.max = R.max,
      RT.opt = RT.opt,
      RT.max = RT.max,  
      M = M,
      a = a,
      c = c,
      v = v
    )
  
  biomass.int.s.consumer.foodweb(tmax, dt, B_init, X_init,  params)
  
}
