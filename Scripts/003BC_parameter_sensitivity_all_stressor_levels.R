#Sensitivity 
library(tidyverse)
library(RColorBrewer)
library(patchwork)
library(ggthemes)

source("Functions/conversion_functions.r")
source("Functions/population_model_analytical.R")
source("Functions/get_stressor_interaction.R")

mytime <- format(Sys.time(), "%Y_%m_%d")

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
time_vect <- seq(0, 250, by = 1)
B_vect <- B_t(params = this_param_set, t = time_vect) 
plot( y = B_vect, x = time_vect)

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


experiment_data2 <- compare_param(params = this_param_set, t = time_vect, T. = templight$T., I. =templight$I.) %>% 
  do.call("rbind", .)

light_for_temp <- 400
temp_for_light <- 40

temperature_data <- rbind(get_stressor_interaction(data = experiment_data2, temp = 36, light = light_for_temp),
                          get_stressor_interaction(data = experiment_data2, temp = 38, light = light_for_temp),
                          get_stressor_interaction(data = experiment_data2, temp = 40, light = light_for_temp),
                          get_stressor_interaction(data = experiment_data2, temp = 42, light = light_for_temp)) %>% within({
                            Days <-  t
                            T. <- factor(T., levels = c(42, 40, 38, 36))
                          })

light_data <- rbind(get_stressor_interaction(data = experiment_data2, temp = temp_for_light, light = 200),
                    get_stressor_interaction(data = experiment_data2, temp = temp_for_light, light = 400),
                    get_stressor_interaction(data = experiment_data2, temp = temp_for_light, light = 600),
                    get_stressor_interaction(data = experiment_data2, temp = temp_for_light, light = 800)) %>% within({
                      Days <- t
                      I. <- factor(I.)
                    })

experiment_data2 <- compare_param(params = this_param_set, t = time_vect, T. = templight$T., I. =templight$I.) %>% 
  do.call("rbind", .)

light_for_temp <- 600
temp_for_light <- 38

temperature_data_1 <- rbind(get_stressor_interaction(data = experiment_data2, temp = 36, light = light_for_temp),
                          get_stressor_interaction(data = experiment_data2, temp = 38, light = light_for_temp),
                          get_stressor_interaction(data = experiment_data2, temp = 40, light = light_for_temp),
                          get_stressor_interaction(data = experiment_data2, temp = 42, light = light_for_temp)) %>% within({
                            Days <-  t
                            T. <- factor(T., levels = c(42, 40, 38, 36))
                          })

light_data_1 <- rbind(get_stressor_interaction(data = experiment_data2, temp = temp_for_light, light = 200),
                    get_stressor_interaction(data = experiment_data2, temp = temp_for_light, light = 400),
                    get_stressor_interaction(data = experiment_data2, temp = temp_for_light, light = 600),
                    get_stressor_interaction(data = experiment_data2, temp = temp_for_light, light = 800)) %>% within({
                      Days <- t
                      I. <- factor(I.)
                    })

experiment_data2 <- compare_param(params = this_param_set, t = time_vect, T. = templight$T., I. =templight$I.) %>% 
  do.call("rbind", .)

light_for_temp <- 800
temp_for_light <- 36

temperature_data_2 <- rbind(get_stressor_interaction(data = experiment_data2, temp = 36, light = light_for_temp),
                            get_stressor_interaction(data = experiment_data2, temp = 38, light = light_for_temp),
                            get_stressor_interaction(data = experiment_data2, temp = 40, light = light_for_temp),
                            get_stressor_interaction(data = experiment_data2, temp = 42, light = light_for_temp)) %>% within({
                              Days <-  t
                              T. <- factor(T., levels = c(42, 40, 38, 36))
                            })

light_data_2 <- rbind(get_stressor_interaction(data = experiment_data2, temp = temp_for_light, light = 200),
                      get_stressor_interaction(data = experiment_data2, temp = temp_for_light, light = 400),
                      get_stressor_interaction(data = experiment_data2, temp = temp_for_light, light = 600),
                      get_stressor_interaction(data = experiment_data2, temp = temp_for_light, light = 800)) %>% within({
                        Days <- t
                        I. <- factor(I.)
                      })



myColors <- brewer.pal(4,"YlOrRd")
names(myColors) <- levels(experiment_data$Stressor)
colScale <- scale_colour_manual(name = "Stressor",values = myColors)

names(myColors) <- rev(levels(temperature_data$T.))
colScale <- scale_colour_manual(name = "Temperature",values = myColors)
K <-  this_param_set$B.max

# LRR plot
s1A <- ggplot(temperature_data, aes(x = Days, y = interact_metric, colour = T.))+
  geom_line(size = 2)+
  colScale+
  ggtitle("Light = 400")+
  geom_hline(yintercept = 0, lty = 2, size = 0.8)+
  xlab("Time (days)") + 
  ylim(-0.5,0.55)+
  theme(plot.title = element_text(hjust = 0.5))+
  ylab(expression(rho)) +
  theme_clean()+
  theme(axis.text.x = element_text(size = 14))+
  theme(axis.text.y = element_text(size = 14))+
  theme(axis.title.x = element_text(size = 18))+
  theme(axis.title.y = element_text(size = 18))
s1A


names(myColors) <- rev(levels(light_data$I.))
colScale <- scale_colour_manual(name = "Light",values = myColors)
K <-  this_param_set$B.max

s1B <- ggplot(light_data, aes(x = Days, y = interact_metric, colour = I.))+
  geom_line(size = 2)+
  colScale+
  ggtitle("Temperature = 40")+
  geom_hline(yintercept = 0, lty = 2, size = 0.8)+
  xlab("Time (days)")  + 
  ylim(-0.5,0.50)+
  theme(plot.title = element_text(hjust = 0.5))+
  ylab(expression(rho)) + labs(color = "Light") +
  theme_clean()+
  theme(axis.text.x = element_text(size = 14))+
  theme(axis.text.y = element_text(size = 14))+
  theme(axis.title.x = element_text(size = 18))+
  theme(axis.title.y = element_text(size = 18))
s1B

# LRR plot
names(myColors) <- rev(levels(temperature_data_1$T.))
colScale <- scale_colour_manual(name = "Temperature",values = myColors)
K <-  this_param_set$B.max
s1C <- ggplot(temperature_data_1, aes(x = Days, y = interact_metric, colour = T.))+
  geom_line(size = 2)+
  colScale+
  ggtitle("Light = 600")+  
  ylim(-0.5,0.55)+
  geom_hline(yintercept = 0, lty = 2, size = 0.8)+
  xlab("Time (days)") + 
  theme(plot.title = element_text(hjust = 0.5))+
  ylab(expression(rho)) +
  theme_clean()+
  theme(axis.text.x = element_text(size = 14))+
  theme(axis.text.y = element_text(size = 14))+
  theme(axis.title.x = element_text(size = 18))+
  theme(axis.title.y = element_text(size = 18))
s1C


names(myColors) <- rev(levels(light_data_1$I.))
colScale <- scale_colour_manual(name = "Light",values = myColors)
K <-  this_param_set$B.max

s1D <- ggplot(light_data_1, aes(x = Days, y = interact_metric, colour = I.))+
  geom_line(size = 2)+
  colScale+
  ggtitle("Temperature = 38")+
  ylim(-0.5,0.5)+
  geom_hline(yintercept = 0, lty = 2, size = 0.8)+
  xlab("Time (days)")  + 
  theme(plot.title = element_text(hjust = 0.5))+
  ylab(expression(rho)) + labs(color = "Light") +
  theme_clean()+
  theme(axis.text.x = element_text(size = 14))+
  theme(axis.text.y = element_text(size = 14))+
  theme(axis.title.x = element_text(size = 18))+
  theme(axis.title.y = element_text(size = 18))
s1D

# LRR plot
names(myColors) <- rev(levels(temperature_data_2$T.))
colScale <- scale_colour_manual(name = "Temperature",values = myColors)
K <-  this_param_set$B.max
s1E <- ggplot(temperature_data_2, aes(x = Days, y = interact_metric, colour = T.))+
  geom_line(size = 2)+
  colScale+
  ggtitle("Light = 800")+
  ylim(-0.5,0.55)+
  geom_hline(yintercept = 0, lty = 2, size = 0.8)+
  xlab("Time (days)") + 
  theme(plot.title = element_text(hjust = 0.5))+
  ylab(expression(rho)) +
  theme_clean()+
  theme(axis.text.x = element_text(size = 14))+
  theme(axis.text.y = element_text(size = 14))+
  theme(axis.title.x = element_text(size = 18))+
  theme(axis.title.y = element_text(size = 18))
s1E


names(myColors) <- rev(levels(light_data_2$I.))
colScale <- scale_colour_manual(name = "Light",values = myColors)
K <-  this_param_set$B.max

s1F <- ggplot(light_data_2, aes(x = Days, y = interact_metric, colour = I.))+
  geom_line(size = 2)+
  colScale+
  ggtitle("Tempperature = 36")+
  ylim(-0.5,0.5)+
  geom_hline(yintercept = 0, lty = 2, size = 0.8)+
  xlab("Time (days)")  + 
  theme(plot.title = element_text(hjust = 0.5))+
  ylab(expression(rho)) + labs(color = "Light") +
  theme_clean()+
  theme(axis.text.x = element_text(size = 14))+
  theme(axis.text.y = element_text(size = 14))+
  theme(axis.title.x = element_text(size = 18))+
  theme(axis.title.y = element_text(size = 18))
s1F

s1 <- (s1A + s1B) / (s1C + s1D) / (s1E +s1F)
s1 <- s1 + plot_annotation(tag_levels = 'A')
s1

ggsave(path = "Plots", filename = paste0(mytime, "_sensitivity.png"), s1, width = 12, height = 8, units = c("in"), dpi = 300)



# CONSUMER - RESOURCE

# setwd and load libraries 
library(tidyverse)
library(patchwork)
library(gridExtra)
library(gganimate)
library(RColorBrewer)
library(ggthemes)

mytime <- format(Sys.time(), "%Y_%m_%d")

# source functions 
source("Functions/conversion_functions.r")
source("Functions/get_stressor_interaction.R")
source("Functions/photo_temp_function.R")
source("Functions/resp_function.R")
source("Functions/photo_light_shade_function.R")
source("Functions/photosyn_shade_function.R")
source("Functions/biomass.integration_shade_function_consumer_foodweb.R")
source("Functions/dBdt_shade_function_consumer_foodweb.R")
source("Functions/attack_mort_temp_scaling_function.R")


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


light_for_temp <- 400
temp_for_light <- 40

temperature_data <- rbind(get_stressor_interaction(data = temp_light_vary, temp = 36, light = light_for_temp),
                          get_stressor_interaction(data = temp_light_vary, temp = 38, light = light_for_temp),
                          get_stressor_interaction(data = temp_light_vary, temp = 40, light = light_for_temp),
                          get_stressor_interaction(data = temp_light_vary, temp = 42, light = light_for_temp)) %>% within({
                            Days <-  (time-1)/24
                            T. <- factor(T., levels = c(42, 40, 38, 36))
                          })

light_data <- rbind(get_stressor_interaction(data = temp_light_vary, temp = temp_for_light, light = 200),
                    get_stressor_interaction(data = temp_light_vary, temp = temp_for_light, light = 400),
                    get_stressor_interaction(data = temp_light_vary, temp = temp_for_light, light = 600),
                    get_stressor_interaction(data = temp_light_vary, temp = temp_for_light, light = 800)) %>% within({
                      Days = (time-1)/24
                      I. <- factor(I., levels = c(200, 400, 600, 800))
                    })

temp_light_vary <- xout3 # match ovs here for light_vary and temp_vary 

light_for_temp <- 600
temp_for_light <- 38

temperature_data_1 <- rbind(get_stressor_interaction(data = temp_light_vary, temp = 36, light = light_for_temp),
                          get_stressor_interaction(data = temp_light_vary, temp = 38, light = light_for_temp),
                          get_stressor_interaction(data = temp_light_vary, temp = 40, light = light_for_temp),
                          get_stressor_interaction(data = temp_light_vary, temp = 42, light = light_for_temp)) %>% within({
                            Days <-  (time-1)/24
                            T. <- factor(T., levels = c(42, 40, 38, 36))
                          })

light_data_1 <- rbind(get_stressor_interaction(data = temp_light_vary, temp = temp_for_light, light = 200),
                    get_stressor_interaction(data = temp_light_vary, temp = temp_for_light, light = 400),
                    get_stressor_interaction(data = temp_light_vary, temp = temp_for_light, light = 600),
                    get_stressor_interaction(data = temp_light_vary, temp = temp_for_light, light = 800)) %>% within({
                      Days = (time-1)/24
                      I. <- factor(I., levels = c(200, 400, 600, 800))
                    })

temp_light_vary <- xout3 # match ovs here for light_vary and temp_vary 

light_for_temp <- 800
temp_for_light <- 36

temperature_data_2 <- rbind(get_stressor_interaction(data = temp_light_vary, temp = 36, light = light_for_temp),
                          get_stressor_interaction(data = temp_light_vary, temp = 38, light = light_for_temp),
                          get_stressor_interaction(data = temp_light_vary, temp = 40, light = light_for_temp),
                          get_stressor_interaction(data = temp_light_vary, temp = 42, light = light_for_temp)) %>% within({
                            Days <-  (time-1)/24
                            T. <- factor(T., levels = c(42, 40, 38, 36))
                          })

light_data_2 <- rbind(get_stressor_interaction(data = temp_light_vary, temp = temp_for_light, light = 200),
                    get_stressor_interaction(data = temp_light_vary, temp = temp_for_light, light = 400),
                    get_stressor_interaction(data = temp_light_vary, temp = temp_for_light, light = 600),
                    get_stressor_interaction(data = temp_light_vary, temp = temp_for_light, light = 800)) %>% within({
                      Days = (time-1)/24
                      I. <- factor(I., levels = c(200, 400, 600, 800))
                    })

names(myColors) <- rev(levels(temperature_data$T.))
colScale <- scale_colour_manual(name = "Temperature",values = myColors)
K <-  this_param_set$B.max

# LRR plot
s2A <- ggplot(temperature_data, aes(x = Days, y = interact_metric, colour = T.))+
  geom_line(size = 2)+
  colScale+
  ggtitle("Light = 400")+
  
  ylim(-.5,.25)+
  geom_hline(yintercept = 0, lty = 2, size = 0.8)+
  xlab("Time (days)") + 
  theme(plot.title = element_text(hjust = 0.5))+
  ylab(expression(rho)) +
  theme_clean()+
  theme(axis.text.x = element_text(size = 14))+
  theme(axis.text.y = element_text(size = 14))+
  theme(axis.title.x = element_text(size = 18))+
  theme(axis.title.y = element_text(size = 18))
s2A


names(myColors) <- rev(levels(light_data$I.))
colScale <- scale_colour_manual(name = "Light",values = myColors)
K <-  this_param_set$B.max

s2B <- ggplot(light_data, aes(x = Days, y = interact_metric, colour = I.))+
  geom_line(size = 2)+
  ggtitle("Temperature = 40")+
  colScale+
  geom_hline(yintercept = 0, lty = 2, size = 0.8)+
  ylim(-.5,.25)+
  xlab("Time (days)")  + 
  theme(plot.title = element_text(hjust = 0.5))+
  ylab(expression(rho)) + labs(color = "Light") +
  theme_clean()+
  theme(axis.text.x = element_text(size = 14))+
  theme(axis.text.y = element_text(size = 14))+
  theme(axis.title.x = element_text(size = 18))+
  theme(axis.title.y = element_text(size = 18))
s2B


names(myColors) <- rev(levels(temperature_data_1$T.))
colScale <- scale_colour_manual(name = "Temperature",values = myColors)
K <-  this_param_set$B.max

# LRR plot
s2C <- ggplot(temperature_data_1, aes(x = Days, y = interact_metric, colour = T.))+
  geom_line(size = 2)+
  colScale+
  ggtitle("Light = 600")+
  ylim(-.5,.25)+
  geom_hline(yintercept = 0, lty = 2, size = 0.8)+
  xlab("Time (days)") + 
  theme(plot.title = element_text(hjust = 0.5))+
  ylab(expression(rho)) +
  theme_clean()+
  theme(axis.text.x = element_text(size = 14))+
  theme(axis.text.y = element_text(size = 14))+
  theme(axis.title.x = element_text(size = 18))+
  theme(axis.title.y = element_text(size = 18))
s2C


names(myColors) <- rev(levels(light_data_1$I.))
colScale <- scale_colour_manual(name = "Light",values = myColors)
K <-  this_param_set$B.max

s2D <- ggplot(light_data_1, aes(x = Days, y = interact_metric, colour = I.))+
  geom_line(size = 2)+
  colScale+
  geom_hline(yintercept = 0, lty = 2, size = 0.8)+
  ylim(-.5,.25)+
  ggtitle("Temperature = 38")+
  xlab("Time (days)")  + 
  theme(plot.title = element_text(hjust = 0.5))+
  ylab(expression(rho)) + labs(color = "Light") +
  theme_clean()+
  theme(axis.text.x = element_text(size = 14))+
  theme(axis.text.y = element_text(size = 14))+
  theme(axis.title.x = element_text(size = 18))+
  theme(axis.title.y = element_text(size = 18))
s2D

names(myColors) <- rev(levels(temperature_data_2$T.))
colScale <- scale_colour_manual(name = "Temperature",values = myColors)
K <-  this_param_set$B.max

# LRR plot
s2E <- ggplot(temperature_data_2, aes(x = Days, y = interact_metric, colour = T.))+
  geom_line(size = 2)+
  colScale+
  ylim(-.5,.25)+
  geom_hline(yintercept = 0, lty = 2, size = 0.8)+
  ggtitle("Light = 800")+
  xlab("Time (days)") + 
  theme(plot.title = element_text(hjust = 0.5))+
  ylab(expression(rho)) +
  theme_clean()+
  theme(axis.text.x = element_text(size = 14))+
  theme(axis.text.y = element_text(size = 14))+
  theme(axis.title.x = element_text(size = 18))+
  theme(axis.title.y = element_text(size = 18))
s2E


names(myColors) <- rev(levels(light_data_2$I.))
colScale <- scale_colour_manual(name = "Light",values = myColors)
K <-  this_param_set$B.max

s2F <- ggplot(light_data_2, aes(x = Days, y = interact_metric, colour = I.))+
  geom_line(size = 2)+
  colScale+
  geom_hline(yintercept = 0, lty = 2, size = 0.8)+
  ggtitle("Temperature = 36")+
  ylim(-.5,.25)+
  xlab("Time (days)")  + 
  theme(plot.title = element_text(hjust = 0.5))+
  ylab(expression(rho)) + labs(color = "Light") +
  theme_clean()+
  theme(axis.text.x = element_text(size = 14))+
  theme(axis.text.y = element_text(size = 14))+
  theme(axis.title.x = element_text(size = 18))+
  theme(axis.title.y = element_text(size = 18))
s2F


s2 <- (s2A + s2B) / (s2C + s2D) / (s2E +s2F)
s2 <- s2 + plot_annotation(tag_levels = 'A')
s2

ggsave(path = "Plots", filename = paste0(mytime, "_parameter_sensitivity_consumner-resource_all-stressors.png"), s2, width = 12, height = 8, units = c("in"), dpi = 300)


