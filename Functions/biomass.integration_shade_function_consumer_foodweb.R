# ---------------------------------
# Biomass integration 
# ---------------------------------
# looks at change in biomass through time based on netgrowth and mortality 
# tmax - length of vector
# dt - timestep
# B1 = initial biomass

# calling dBdt.s function with calculates change in biomass based on net growth - mortality 
biomass.int.s.consumer.foodweb <- function(params, model){

          B <- rep(NA, params$tmax) # create empyt vector to store output
          B[1] <- params$B_init         # set inital biomass in time 1
  
          X <- rep(NA, params$tmax) # create empty vector to store output
          X[1] <- params$X_init
          
          # retrieve the model of interest for solving numerically
          dBdt_mod <- switch(model, "base" = dBdt.s.consumer.foodweb,
                             "gompertz" = dBdt.gompertz,
                             "pmaxdepend" = dBdt.pmaxdepend,
                             "holling2" = dBdt.holling2)
  
          if (length(dBdt_mod) == 0) { stop("Model choice is invalid - should be one of:
          base, gompertz, pmaxdepend, holling2")}
            
  for (time_step in 1:(params$tmax-1)){
        dBdX <- dBdt_mod(1, B[time_step], X[time_step], params)
          B[time_step+1] <- B[time_step] + dBdX[1]*params$dt
          X[time_step+1] <- X[time_step] + dBdX[2]*params$dt
      }
 mylist <- list(B,X)
 return(mylist)
}

# create wrapper to list parameters before calling pmap 
biomass.int.wrapper.cons.foodweb <- function(..., model = "base"){ 
  
  params <- list(...)
  biomass.int.s.consumer.foodweb(params, model)
  
}
