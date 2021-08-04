
# ------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------
#
# Phase plots for BIOMASS FOODWEB MODEL 
#
# ------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------

source(file = "Functions/MC_conversion_functions.R")
source(file = "Functions/resp_function.R")
source(file = "Functions/MC_attack_mort_temp_scaling_function.R")

library(phaseR)


# Input required parameters  ----------------------------------------------

this_param_set <- list(
  B.max = 667, 
  T. = 40,      # Temp                          degrees C
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
  v = 0.05     # consumer mortality 
)



# Generate the phase plot -------------------------------------------------


# written in the format for phaseR package
dBdt.s.consumer.foodweb <- function(t, y, parameters){
  with(as.list(parameters), { 
    
    B <- y[1]
    X <- y[2]
    
    phi <-  tanh(I./Ik)*PT.max*(T.max-T.)/(T.max - T.opt) * 
      (T./T.opt)^(T.opt/(T.max - T.opt))
    R <- Resp.temp(R.max, T., RT.opt, RT.max)
    S <- Attack.Mortal.Temp.Scaling(T.)
    
    dy <- numeric(2)
    dy[1] <- (phi*(1-B/B.max) - R - M - (a*S*X))*B
    dy[2] <- (a*S*c*B*X)-(v*S*X)
    list(dy)
  })
}

# What stressors to we want to compare
tempsens <- c(35, 42) 
lightsens <- c(1000, 200)
templight <- expand.grid(T. = tempsens, I. = lightsens)

treatment_name <- c("Control", "Light Stress", "Temp Stress", 
                    "Temp + Light Stress")

# Specify matrix of initial conditions
#y0 <- matrix(c(66, 1, 500, 20, 200, 35, 500, 80, 300, 100), ncol = 2, nrow = 5, byrow = TRUE) # base set
y0 <- as.matrix(expand.grid(seq(1,580, length.out = 10), seq(1,130, length.out = 10)))

png(filename = "Plots/phase_plots4.png", width = 22, height = 22, units = "cm",  res = 300)
par(mfrow = c(2,2))
par(mar = c(4.5, 4.8, 1.1, 0.9))

for (i in 1:nrow(templight)){
  
  # Update the parameter set
  this_param_set[c("T.", "I.")] <- templight[i,]

  { # Generate the phase plot for the given parameters
  lotkaVolterra.flowField <- flowField(dBdt.s.consumer.foodweb, 
                                       xlim = c(0, 600), ylim = c(0, 130),
                                       state.names = c("B","X"),
                                       parameters = this_param_set, points = 30,
                                       add = FALSE, cex.lab = 1.5, yaxt = "none", xaxt = "none", 
                                       xlab = bquote('Biomass ('*'g' ~DW ~m^-2*')'), 
                                       ylab = expression(Consumer ~density ~(m^-2)))
  
  lotkaVolterra.nullclines <- nullclines(dBdt.s.consumer.foodweb, 
                                         xlim = c(-1, 600), ylim = c(-1, 130),
                                         parameters = this_param_set, points = 500,
                                         col = c("blue", "red"),
                                         add.legend = FALSE, cex = 2.5)
  
  lotkaVolterra.trajectory <- trajectory(dBdt.s.consumer.foodweb, y0 = y0, tlim = c(0,1000), tstep = 1/24,
                                         parameters = this_param_set, colour = rep("black", nrow(y0)),
                                         state.names = c("B","X"))
  
  text(x = 25,y = 120, LETTERS[[i]], cex = 2)
  axis(1, at = seq(0,600, by = 100), las=1, cex.axis = 1.3)
  axis(2, at = seq(0,120, by = 20), las=1, cex.axis = 1.3)
  #title(xlab=bquote('Biomass ('*'g' ~DW ~m^-2*')'), ylab="y-axis label") 
  
  #legend(y = 120, x = 450, legend = c("B nullcline", "X nullcline"), 
         #lty=c(1,1), col = c("blue", "red"), cex=0.5 )
  #title(main = treatment_name[[i]])
  }

}
dev.off()

