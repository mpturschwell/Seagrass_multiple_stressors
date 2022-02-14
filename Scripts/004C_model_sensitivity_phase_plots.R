
# ------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------
#
# Phase plots for BIOMASS FOODWEB MODEL 
#
# ------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------

source(file = "Functions/conversion_functions.R")
source(file = "Functions/resp_function.R")
source(file = "Functions/attack_mort_temp_scaling_function.R")

library(phaseR)
library(magrittr)


# Input required parameters  ----------------------------------------------

this_param_set <- list(
  B.max = 667, 
  T. = 35,      # Temp                          degrees C
  I. = 1000,     # Irradiance                    mmol m-2 second-1
  PT.max = mg.C.per.h_to_g.dry.wt.per.d(4.7),   # Max gross production (temp)   mg C g-1 dry wt h-1
  PL.max = umol_to_mg(422), # Max gross production (light)  mmol O2 g dry wt-1 hr-1
  T.opt = 34.9, # optimum temperature           degrees C
  T.max = 44.5, # Maximum temperature           degrees C
  Ik = 319,     # Saturation irradiance         mmol m-2 second-1
  R.max = mg.C.per.h_to_g.dry.wt.per.d(1.1),  # Maximum respiration           mg C g-1 dry wt h-1
  RT.opt = 39.1,  # Resp optimum temperature      degrees C
  RT.max = 45.6,
  M = 0.004,
  tmax = 730,   # length of simulation 
  dt = 1,     # time interval
  B_init = 66,    # Inital Seagrass proportional biomass
  a = 0.001,     # attack rate 
  X_init = 1,    # consumer initial density 
  c = 0.15,    # Assimilation efficiency (10-20%)
  v = 0.05,#0.05,     # consumer mortality 
  h = 3      # Handling time for the holling type II model
  
) %>% within({
  # Parameter for the Pmax dependent model
  psi <-  tanh(I./Ik)*PT.max*(T.max-T.)/(T.max - T.opt) * 
               (T./T.opt)^(T.opt/(T.max - T.opt))/B.max
  })



# Generate the phase plot -------------------------------------------------


# written in the format for phaseR package
dBdt_base_mod<- function(t, y, parameters){
  with(as.list(parameters), { 
    
    B <- y[1]
    X <- y[2]
    
    P.max <-  tanh(I./Ik)*PT.max*(T.max-T.)/(T.max - T.opt) * 
      (T./T.opt)^(T.opt/(T.max - T.opt))
    R <- Resp.temp(R.max, T., RT.opt, RT.max)
    S <- Attack.Mortal.Temp.Scaling(T.)
    
    dy <- numeric(2)
    dy[1] <- (P.max*(1-B/B.max) - R - M - (a*S*X))*B
    dy[2] <- (a*S*c*B*X)-(v*S*X)
    list(dy)
  })
}

dBdt_gompertz_mod<- function(t, y, parameters){
  with(as.list(parameters), { 
    
    B <- y[1]
    X <- y[2]
    
    P.max <-  tanh(I./Ik)*PT.max*(T.max-T.)/(T.max - T.opt) * 
      (T./T.opt)^(T.opt/(T.max - T.opt))
    R <- Resp.temp(R.max, T., RT.opt, RT.max)
    S <- Attack.Mortal.Temp.Scaling(T.)
    
    dy <- numeric(2)
    dy[1] <- (P.max*log(B.max/(B+1e-8)) - R - M - (a*S*X))*B # This is a hack to make it plot
    dy[2] <- (a*S*c*B*X)-(v*S*X)
    list(dy)
  })
}

dBdt_pmaxdep_mod<- function(t, y, parameters){
  with(as.list(parameters), { 
    
    B <- y[1]
    X <- y[2]
    
    P.max <-  tanh(I./Ik)*PT.max*(T.max-T.)/(T.max - T.opt) * 
      (T./T.opt)^(T.opt/(T.max - T.opt))
    R <- Resp.temp(R.max, T., RT.opt, RT.max)
    S <- Attack.Mortal.Temp.Scaling(T.)
    
    dy <- numeric(2)
    dy[1] <- ((P.max-psi*B) - R - M - (a*S*X))*B
    dy[2] <- (a*S*c*B*X)-(v*S*X)
    list(dy)
  })
}

dBdt_holling2_mod<- function(t, y, parameters){
  with(as.list(parameters), { 
    
    B <- y[1]
    X <- y[2]
    
    P.max <-  tanh(I./Ik)*PT.max*(T.max-T.)/(T.max - T.opt) * 
      (T./T.opt)^(T.opt/(T.max - T.opt))
    R <- Resp.temp(R.max, T., RT.opt, RT.max)
    S <- Attack.Mortal.Temp.Scaling(T.)
    
    dy <- numeric(2)
    dy[1] <- (P.max*(1-B/B.max) - R - M)*B - (a*S*B)/(1 + (a*h*S*B))*X
    dy[2] <- (a*S*B)/(1 + (a*h*S*B))*c*X - v*S*X
    list(dy)
  })
}

model_vect <- list(dBdt_base_mod, dBdt_gompertz_mod, dBdt_pmaxdep_mod, dBdt_holling2_mod)

# Specify matrix of initial conditions
#y0 <- matrix(c(66, 1, 500, 20, 200, 35, 500, 80, 300, 100), ncol = 2, nrow = 5, byrow = TRUE) # base set
y0 <- as.matrix(expand.grid(seq(40,560, length.out = 4), seq(20,280, length.out = 4)))

# Specify stressor scenarios
templight <- expand.grid(T. = c(35, 42) , I. = c(1000, 200))

for (j in seq_along(model_vect)){
  
  file_name <- sprintf("Plots/%s_phase_plots.png", c("base_model", "gompertz", 
                                                     "pmaxdepend", "hollingII")[[j]])
  mod <- model_vect[[j]]
  
  png(filename = file_name, width = 22, height = 22, units = "cm",  res = 300)
  par(mfrow = c(2,2))
  par(mar = c(4.5, 4.8, 1.1, 0.9))

for (i in seq_along(model_vect)){
  
  # Make the parameter set of interest
  temp_params <- within(this_param_set,{
    
    # Change mortality rate of consumer due to handling time in holling
    if (j == 4) v <- 0.0275 
    
    # Specify stressor scenario
    T. <- templight$T.[[i]]
    I. <- templight$I.[[i]]
    
    })
  
                        
  { # Generate the phase plot for the given parameters
  lotkaVolterra.flowField <- flowField(mod, 
                                       xlim = c(0, 600), ylim = c(0, 300),
                                       state.names = c("B","X"),
                                       parameters = temp_params, points = 30,
                                       add = FALSE, cex.lab = 1.5, yaxt = "none", xaxt = "none", 
                                       xlab = bquote('Biomass ('*'g' ~DW ~m^-2*')'), 
                                       ylab = expression(Consumer ~density ~(m^-2)))
  
  lotkaVolterra.nullclines <- nullclines(mod, 
                                         xlim = c(-1, 600), ylim = c(-1, 300),
                                         parameters = temp_params, points = 500,
                                         col = c("blue", "red"),
                                         add.legend = FALSE, cex = 2.5)
  
  lotkaVolterra.trajectory <- trajectory(mod, y0 = y0, tlim = c(0,1000), tstep = 1/24,
                                         parameters = temp_params, colour = rep("black", nrow(y0)),
                                         state.names = c("B","X"))

  text(x = 25,y = 300, LETTERS[[i]], cex = 2)
  axis(1, at = seq(0,600, by = 100), las=1, cex.axis = 1.3)
  axis(2, at = seq(0,300, by = 50), las=1, cex.axis = 1.3)
  #title(xlab=bquote('Biomass ('*'g' ~DW ~m^-2*')'), ylab="y-axis label") 
  
  #legend(y = 120, x = 450, legend = c("B nullcline", "X nullcline"), 
         #lty=c(1,1), col = c("blue", "red"), cex=0.5 )
  #title(main = treatment_name[[i]])
  }

}
dev.off()
}
