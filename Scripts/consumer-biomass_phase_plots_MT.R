
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
#only select control and multiple stressor 
templight <- templight[c(1,4),]
treatment_name <- c("Control", "Light Stress", "Temp Stress", 
                    "Temp + Light Stress")

# Specify matrix of initial conditions
#y0 <- matrix(c(66, 1, 500, 20, 200, 35, 500, 80, 300, 100), ncol = 2, nrow = 5, byrow = TRUE) # base set
y0 <- as.matrix(expand.grid(seq(1,580, length.out = 5), seq(1,130, length.out = 5)))

png(filename = "Plots/phase_plots4_mt_AB.png", width = 22, height = 11, units = "cm",  res = 300)
par(mfrow = c(1,2))
par(mar = c(4.5, 4.8, 1.1, 0.9))

for (i in 1:nrow(templight)){
  # Update the parameter set
#fig1A  
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
                                         col = c("blue", "red"), lwd = 3,
                                         add.legend = FALSE, cex = 2.5)
  
  lotkaVolterra.trajectory <- trajectory(dBdt.s.consumer.foodweb, y0 = y0, tlim = c(0,1000), tstep = 1/24,
                                         parameters = this_param_set, colour = rep("black", nrow(y0)),
                                         state.names = c("B","X"))
  
  text(x = 25,y = 120, LETTERS[[i]], cex = 2)
  axis(1, at = seq(0,600, by = 100), las=1, cex.axis = 1.3)
  axis(2, at = seq(0,120, by = 20), las=1, cex.axis = 1.3)

  }
}
dev.off()

# ------------------------------------------------------------------------------------------------------
#   Sensitivity for BIOMASS FOODWEB MODEL 
# ------------------------------------------------------------------------------------------------------

# setwd and load libraries 
library(purrr)
library(ggplot2)
library(tidyr)
library(patchwork)
library(dplyr)
library(gridExtra)
library(gganimate)
library(RColorBrewer)
library(ggthemes)

mytime <- format(Sys.time(), "%Y_%m_%d")

# source functions 
source("Functions/MC_conversion_functions.r")
source("Functions/get_stressor_interaction.R")
source("Functions/photo_temp_function.R")
source("Functions/resp_function.R")
source("Functions/photo_light_shade_function.R")
source("Functions/photosyn_shade_function.R")
#source("Functions/net.growth_shade_function.R")
source("Functions/Foodweb-functions/biomass.integration_shade_function_consumer_foodweb.R")
source("Functions/Foodweb-functions/MC_dBdt_shade_function_consumer_foodweb.R")
source("Functions/MC_attack_mort_temp_scaling_function.R")


# -------------------------------------
# Input required parameters 
# -------------------------------------
this_param_set <- list(
  B.max = 667, 
  T. = 40,      # Temp                          degrees C
  I. = 1000,     # Irradiance                    mmol m-2 second-1
  PT.max = mg.C.per.h_to_g.dry.wt.per.d(4.7),  # Max gross production (temp)   mg C g-1 dry wt h-1
  PL.max = umol_to_mg(422), # Max gross production (light)  mmol O2 g dry wt-1 hr-1
  T.opt = 34.9, # optimum temperature           degrees C
  T.max = 44.5, # Maximum temperature           degrees C
  Ik = 319,     # Saturation irradiance         mmol m-2 second-1
  R.max = mg.C.per.h_to_g.dry.wt.per.d(1.1), # Maximum respiration           mg C g-1 dry wt h-1
  RT.opt = 39.1,  # Resp optimum temperature      degrees C
  RT.max = 45.6,
  M = 0.004, #0.004,
  tmax = 2*365*24,   # length of simulation 
  dt = 1/24,     # time interval
  B_init = 66,    # Inital Seagrass proportional biomass
  a = 0.001,     # attack rate 0.001 
  X_init = 1,    # consumer initial density 
  c = 0.15,    # Assimilation efficiency (10-20%)
  v = 0.05     # consumer mortality 
)


# testing under optimal conditions 
xout <- pmap(this_param_set, biomass.int.wrapper.cons.foodweb) 
plot(y= xout[1][[1]][[2]], 
     x = seq(this_param_set[["dt"]], 
             this_param_set[["tmax"]]*this_param_set[["dt"]], by = this_param_set[["dt"]]), 
     type = "l", col = "red", ylim = c(0,667),lwd = 2, xlab = "Days", 
     ylab = "Biomass (grams of seagrass tissue per day (m2)")
lines(y = xout[1][[1]][[1]], 
      x = seq(this_param_set[["dt"]], this_param_set[["tmax"]]*this_param_set[["dt"]], 
              by = this_param_set[["dt"]]), col = "blue", lwd = 2)

# ---------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------
# create parameter sequence to run model on and use expand.grid to use all combinations 
# ---------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------


# Generate the treatments we want to simulate
tempsens <- c(35, seq(36, 42, by = 2))
lightsens <- seq(200, 1000, by = 200)
Bsens <- c(145.75, 435.25)
Xsens <- c(33.25, 97.75)
templight <- expand.grid(T. = tempsens, I. = lightsens, B_init = Bsens, X_init = Xsens)

# Remove the treatments we don't need
#templight <- templight[templight$T. %in% c(35, 42.58) | templight$I. %in% c(254.55, 1000),] 

this_param_set$T. <- templight$T.
this_param_set$I. <- templight$I.
this_param_set$B_init <- templight$B_init
this_param_set$X_init <- templight$X_init


xout <- pmap(this_param_set, biomass.int.wrapper.cons.foodweb) 

# extract seagrass biomass vectors only 
xout2 <- lapply(xout, function(x) x[1][[1]])

xout3 <- xout2 %>%
  do.call("rbind", .) %>%
  data.frame() %>%
  cbind(templight) %>%
  tidyr::gather(time, biomass, -T., -I., -B_init, -X_init) %>%
  dplyr::mutate(time = as.numeric(substring(time, 2)))

#templight$net.growth <- xout
temp_light_vary <- xout3 # match ovs here for light_vary and temp_vary 


# Stressor interaction plots ----------------------------------------------
fwm <- temp_light_vary
temp_light_vary <- fwm %>% filter(B_init == 145.75, X_init == 33.25)
temp_light_vary1 <- fwm %>% filter(B_init == 435.25, X_init == 33.25)
temp_light_vary2 <- fwm %>% filter(B_init == 145.75, X_init == 97.75)
temp_light_vary3 <- fwm %>% filter(B_init == 435.25, X_init == 97.75)

init <- get_stressor_interaction(data= temp_light_vary, temp = 42, light = 200) %>% mutate(init = '2')
init1 <- get_stressor_interaction(data= temp_light_vary1, temp = 42, light = 200)%>% mutate(init = '1')
init2 <- get_stressor_interaction(data= temp_light_vary2, temp = 42, light = 200)%>% mutate(init = '4')
init3 <- get_stressor_interaction(data= temp_light_vary3, temp = 42, light = 200)%>% mutate(init = '3')

out <- rbind(init, init1, init2, init3) %>% within({
  Days = (time-1)/24})
out$B_init <- factor(out$init)

K <-  this_param_set$B.max


library(RColorBrewer)
myColors <- brewer.pal(4,"Dark2")
names(myColors) <- levels(out$init)
colScale <- scale_colour_manual(name = "Initial Conditions",values = myColors)

# LRR plot
g1 <- ggplot(out, aes(x = Days, y = interact_metric, colour = init))+
  geom_line(size = 2)+
colScale+
  ylim(-1.85,0.25)+
  geom_hline(yintercept = 0, lty = 2, size = 0.8)+
  xlab("Time (days)") + 
  theme(plot.title = element_text(hjust = 0.5))+
  ylab(expression(rho)) +
  theme_clean()+
  theme(axis.text.x = element_text(size = 14))+
  theme(axis.text.y = element_text(size = 14))+
  theme(axis.title.x = element_text(size = 18))+
  theme(axis.title.y = element_text(size = 18))+
  annotate("text", x=-40, y=0.23, label= "C",  size = 10)  
g1


ggsave(path = "Plots", filename = paste0(mytime, "_Figure5_FOODWEB_TEMP-LIGHT_initial_conditions_1.png"), g1, width = 22, height = 11, units = c("cm"), dpi = 300)


