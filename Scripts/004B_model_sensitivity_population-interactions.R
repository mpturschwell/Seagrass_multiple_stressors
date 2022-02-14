
library(tidyverse)
library(RColorBrewer)
library(patchwork)

source("Functions/conversion_functions.R")
source("Functions/population_model_analytical.R")
source("Functions/get_stressor_interaction.R")

# Function for compare treatments for a given biomass function and params
compare_param <- function(params, t, B_t_fun = B_t, ...){
  # ================================
  # Specify any number of parameter changes
  # wish to see sensitivity after the params and t arguments
  # ================================
  
  # Make all combinations of parameter changes we wish to test
  change_pars <- expand.grid(list(...))
  
  # Find the analytical solution for a single row of parameter changes
  B_t_row <- function(x){
    params[names(change_pars)] <- x # Change params
    row_data <- data.frame(t = t, B = B_t_fun(params = params, t = t)) # Compute solutions for timesteps
    
    # Join on parameters
    for (i in seq_along(names(change_pars))) row_data[[names(change_pars)[[i]]]] <- x[[i]] 
    
    return(row_data)
  }
  
  # Find all rows of parameter changes
  apply(change_pars, MARGIN = 1, FUN = B_t_row )  
  
}  

# Make an expression of the code we want to re-evaluate for different models 
# ===========================================================================
make_plots_expr <- expr({
  
  light_for_temp <- 200
  temp_for_light <- 42
  
  temperature_data <- rbind(get_stressor_interaction(data = experiment_data2, temp = 20, light = light_for_temp),
                            get_stressor_interaction(data = experiment_data2, temp = 30, light = light_for_temp),
                            get_stressor_interaction(data = experiment_data2, temp = 40, light = light_for_temp),
                            get_stressor_interaction(data = experiment_data2, temp = 42, light = light_for_temp)) %>% within({
                              Days <-  t
                              T. <- factor(T., levels = c(42, 40, 30, 20))
                            })
  
  light_data <- rbind(get_stressor_interaction(data = experiment_data2, temp = temp_for_light, light = 200),
                      get_stressor_interaction(data = experiment_data2, temp = temp_for_light, light = 400),
                      get_stressor_interaction(data = experiment_data2, temp = temp_for_light, light = 600),
                      get_stressor_interaction(data = experiment_data2, temp = temp_for_light, light = 800)) %>% within({
                        Days <- t
                        I. <- factor(I.)
                      })
  
  myColors <- brewer.pal(4,"YlOrRd")
  names(myColors) <- rev(levels(temperature_data$T.))
  colScale <- scale_colour_manual(name = "Temperature",values = myColors)
  K <-  this_param_set$B.max
  
  # LRR plot
  g4A <- ggplot(temperature_data, aes(x = Days, y = interact_metric, colour = T.))+
    geom_line(size = 2)+
    colScale+
    #  ylim(-1.25,1.25)+
    geom_hline(yintercept = 0, lty = 2, size = 0.8)+
    xlab("Time (days)") + 
    theme(plot.title = element_text(hjust = 0.5))+
    ylab(expression(I[R])) +
    theme_bw()
  g4A
  
  
  names(myColors) <- rev(levels(light_data$I.))
  colScale <- scale_colour_manual(name = "Light",values = myColors)
  K <-  this_param_set$B.max
  
  g4B <- ggplot(light_data, aes(x = Days, y = interact_metric, colour = I.))+
    geom_line(size = 2)+
    colScale+
    geom_hline(yintercept = 0, lty = 2, size = 0.8)+
    # ylim(-1.25,1.25)+
    xlab("Time (days)")  + 
    theme(plot.title = element_text(hjust = 0.5))+
    ylab(expression(I[R])) + labs(color = "Light") +
    theme_bw()

})

# Specify base parameter set ------------------------------------------------------

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
  B_init = 66    # Initial Seagrass proportional biomass
)

# Generate the treatments we want to compute
tempsens <- c(18:42)
lightsens <- seq(200, 1000, by = 200)
templight <- expand.grid(T. = tempsens, I. = lightsens)

# CHOOSE times
time_vect <- seq(0, 365, by = 1)

# Base model plots --------------------------------------------------------

experiment_data2 <- compare_param(params = this_param_set, t = time_vect,
                                  B_t_fun = B_t, # Specify function here
                                  T. = templight$T., I. =templight$I.) %>% 
  do.call("rbind", .)

eval(make_plots_expr) # Evaluate the expression to generate plots

# bind plots to another name for when we recompute the expression
base_temp_plot <- g4A
base_light_plot <- g4B


# Gompertz model plots ----------------------------------------------------

experiment_data2 <- compare_param(params = this_param_set, t = time_vect,
                                  B_t_fun = B_t_gompertz, # Specify function here
                                  T. = templight$T., I. =templight$I.) %>% 
  do.call("rbind", .)

eval(make_plots_expr) # Evaluate the expression to generate plots

# bind plots to another name for when we recompute the expression
gomp_temp_plot <- g4A
gomp_light_plot <- g4B

# Pmax dependent model plots ----------------------------------------------

source("Functions/photo_temp_function.R")
# Changes to parameter set
pdep_params <- within(this_param_set,{
  psi <- tanh(1000/Ik)*Photo.Temp(T.= 35, PT.max, T.opt, T.max)/B.max
  })

experiment_data2 <- compare_param(params = pdep_params, t = time_vect,
                                  B_t_fun = B_t_pmaxdepend, # Specify function here
                                  T. = templight$T., I. =templight$I.) %>% 
  do.call("rbind", .)
eval(make_plots_expr) # Evaluate the expression to generate plots

# bind plots to another name for when we recompute the expression
pdep_temp_plot <- g4A
pdep_light_plot <- g4B

# Function to change the time limit to the first 50 days
ltrans <- function(x){
  x + xlim(c(0, 50))
  }

# Summary of model plots --------------------------------------------------

# Compare temperature findings
(base_temp_plot + gomp_temp_plot + pdep_temp_plot)/
(base_light_plot + gomp_light_plot + pdep_light_plot)

ggsave(filename = "Plots/model_sens_population_models.png", width = 30, height = 15, units = "cm")

# Compare temperature findings
(ltrans(base_temp_plot) + ltrans(gomp_temp_plot) + ltrans(pdep_temp_plot))/
  (ltrans(base_light_plot) + ltrans(gomp_light_plot) + ltrans(pdep_light_plot))

ggsave(filename = "Plots/transient_model_sens_population_models.png", width = 30, height = 15, units = "cm")
