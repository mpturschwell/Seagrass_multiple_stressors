#fig 3

source("Functions/MC_conversion_functions.r")
source("Functions/population_model_analytical.R")
source("Functions/get_stressor_interaction.R")

mytime <- format(Sys.time(), "%Y_%m_%d")

# Specify base parameter set ------------------------------------------------------

this_param_set <- list(
  B.max = 667, 
  T. = 30,      # Temp                          degrees C
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
  B_init = 66    # Initial Seagrass proportional biomass
)


# Compute the biomass -----------------------------------------------------

# CHOOSE times
#time_vect <- seq(0, 365, by = 1)

time_vect <- seq(0, 250, by = 1)
# Compute B
B_vect <- B_t(params = this_param_set, t = time_vect) 

plot( y = B_vect, x = time_vect)


# Compare treatments ------------------------------------------------------

library(tidyverse)
library(RColorBrewer)
library(patchwork)
library(ggthemes)

# Generate the treatments we want to simulate
tempsens <- c(18:42)
lightsens <- seq(200, 1000, by = 200)
templight <- expand.grid(T. = tempsens, I. = lightsens)


compare_param <- function(params, t, ...){
  # ================================
  # Specify any number of parameter changes
  # wish to see sensitivity after the params and t arguments
  # ================================
  
  # Make all combinations of parameter changes we wish to test
  change_pars <- expand.grid(list(...))
  
  # Find the analytical solution for a single row of parameter changes
  B_t_row <- function(x){
    params[names(change_pars)] <- x # Change params
    row_data <- data.frame(t = t, B = B_t(params = params, t = t)) # Compute solutions for timesteps
    
    # Join on parameters
    for (i in seq_along(names(change_pars))) row_data[[names(change_pars)[[i]]]] <- x[[i]] 
    
    return(row_data)
  }
  
  # Find all rows of parameter changes
  apply(change_pars, MARGIN = 1, FUN = B_t_row )  
  
}  

experiment_data <- compare_param(params = this_param_set, t = time_vect, T. = c(35, 42), I. =c(1000, 200)) %>% 
  do.call("rbind", .) %>% mutate(Stressor = ordered((case_when(T. == 35 & I. == 1000 ~ "Control",
                                                               T. == 35 & I. == 200 ~ "Light Stress",
                                                               T. == 42 & I. == 1000 ~ "Temp Stress",
                                                               T. == 42 & I. == 200 ~ "Temp + Light Stress")), 
                                                    levels = c("Control", "Light Stress", "Temp Stress", 
                                                               "Temp + Light Stress")))


# Mischa's biomass plots ----------------------------------------------------------

#myColors <- brewer.pal(4,"YlOrRd")
myColours <- c("#020006", "#388416", "#095fad", "#E31A1C")
names(myColours) <- levels(experiment_data$Stressor)
colScale <- scale_colour_manual(name = "Stressor",values = myColours)

B.max <-  this_param_set$B.max
g1B <- ggplot(experiment_data, aes(x = t, y = B, colour = Stressor))+
  geom_line(size = 2)+
  colScale+
  xlab("Time (days)") + #labs(colour = "Stressor") +
  theme(plot.title = element_text(hjust = 0.5))+
  ylab(bquote('Biomass ('*'g' ~DW ~m^-2*')')) +
  theme_bw()+ scale_y_continuous(breaks = seq(100, 600, 100)) + 
  geom_hline(yintercept=B.max, lty = 2, size = 1.5)+
  theme_clean()+
  ggtitle("Population sub-model")+
  theme(axis.text.x = element_text(size = 14))+
  theme(axis.text.y = element_text(size = 14))+
  theme(axis.title.x = element_text(size = 14))+
  theme(axis.title.y = element_text(size = 14))

g1B



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
templight <- expand.grid(T. = tempsens, I. = lightsens)

# Remove the treatments we don't need
#templight <- templight[templight$T. %in% c(35, 42.58) | templight$I. %in% c(254.55, 1000),] 

this_param_set$T. <- templight$T.
this_param_set$I. <- templight$I.


xout <- pmap(this_param_set, biomass.int.wrapper.cons.foodweb) 

# extract seagrass biomass vectors only 
xout2 <- lapply(xout, function(x) x[1][[1]])

xout3 <- xout2 %>%
  do.call("rbind", .) %>%
  data.frame() %>%
  cbind(templight) %>%
  tidyr::gather(time, biomass, -T., -I.) %>%
  dplyr::mutate(time = as.numeric(substring(time, 2)))

#templight$net.growth <- xout
temp_light_vary <- xout3 # match ovs here for light_vary and temp_vary 


# ----------------------------------------
# Plot outputs and investigation of metrics
# ---------------------------------------

# make labels 
control <- temp_light_vary %>% 
  filter(T. == 35, I. == 1000) %>% 
  mutate(Stressor = "Control")

Light <- temp_light_vary %>% 
  filter(T. == 35, I. == 200) %>% 
  mutate(Stressor = "Light Stress")

Temp <- temp_light_vary %>% 
  filter(T. == 42, I. == 1000) %>% 
  mutate(Stressor = "Temp Stress")

stress <- temp_light_vary %>% 
  filter(T. == 42, I. == 200) %>% 
  mutate(Stressor = "Temp + Light Stress")

Out <- rbind(control, Light, Temp, stress) %>% within({
  T.<- as.factor(T.)
  I.<- as.factor(I.)
  Stressor <- ordered(Stressor, levels = c("Control", "Light Stress", "Temp Stress", "Temp + Light Stress"))
})

# # Investigation of IR
# IR_data<- data.frame(
#   IR = log((stress$biomass - control$biomass)/(Light$biomass + Temp$biomass - 2*control$biomass)),
#   base_additive = (stress$biomass - control$biomass) - (Light$biomass + Temp$biomass - 2*control$biomass),
#   log_base_additive = (log(stress$biomass) - log(control$biomass)) - (log(Light$biomass) + log(Temp$biomass) - 2*log(control$biomass)),
#   Days = (stress$time-1)/24
# )
# 
# par(mfrow = c(1,2))
# plot(y = IR_data$base_additive, x = IR_data$Days, xlab = "Days", ylab = expression(h[AB]-(h[A]+h[B])))
# plot(y = IR_data$log_base_additive, x = IR_data$Days, xlab = "Days", ylab = expression(log~version~h[AB]-(h[A]+h[B])))
# plot(y = IR_data$IR[-1], x = IR_data$Days[-1], xlab = "Days", ylab = expression(I[R]))


# Biomass over time -------------------------------------------------------


library(RColorBrewer)
#myColors <- brewer.pal(4,"YlOrRd")
myColours <- c("#020006", "#388416", "#095fad", "#E31A1C")
names(myColours) <- levels(Out$Stressor)
colScale <- scale_colour_manual(name = "Stressor",values = myColours)

K <-  this_param_set$B.max
g1C <- ggplot(Out, aes(x = time*this_param_set[["dt"]], y = biomass, colour = Stressor))+
  geom_line(size = 2)+
  colScale+
  xlab("Time (days)") + #labs(colour = "Stressor") + 
  theme(plot.title = element_text(hjust = 0.5))+
  ylab(bquote('Biomass ('*'g' ~DW ~m^-2*')')) +
  theme_bw()+
  geom_hline(yintercept=K, lty = 2, size = 1.5)+
  scale_y_continuous(breaks = seq(100, 600, 100))+
  theme_clean()+
  ggtitle("Consumer-resource model")+
  theme(axis.text.x = element_text(size = 14))+
  theme(axis.text.y = element_text(size = 14))+
  theme(axis.title.x = element_text(size = 14))+
  theme(axis.title.y = element_text(size = 14))

g1C

Figure3 <- g1B/g1C
Figure3 <- Figure2 + plot_annotation(tag_levels = 'A')
Figure3
ggsave(path = "Plots", filename = paste0(mytime, "_Figure3.tiff"), Figure3 , width = 8, height = 8, units = c("in"), dpi = 600)

