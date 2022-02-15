

# setwd and load libraries 
# setwd("/Users/s2904436/Dropbox/Multi-stressor/Code")
library(tidyverse)
library(patchwork)
library(gridExtra)
library(gganimate)
library(RColorBrewer)


# source functions 
source("Functions/conversion_functions.R")
source("Functions/get_stressor_interaction.R")
source("Functions/photo_temp_function.R")
source("Functions/resp_function.R")
source("Functions/photo_light_shade_function.R")
source("Functions/photosyn_shade_function.R")
source("Functions/biomass.integration_shade_function_consumer_foodweb.R")
source("Functions/dBdt_shade_function_consumer_foodweb.R")
source("Functions/alternative_consumer_models_dBdt.R") # alternative models
source("Functions/attack_mort_temp_scaling_function.R")


# ---------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------
# Make an expression of the code we want to re-evaluate for different models 
# ---------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------


# We use an expression so we do not have to keep repeating the code
make_plots_expr <- expr({
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

# Stressor interaction plots ----------------------------------------------


light_for_temp <- 200
temp_for_light <- 42

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


myColors <- brewer.pal(4,"YlOrRd")
names(myColors) <- rev(levels(temperature_data$T.))
colScale <- scale_colour_manual(name = "Temperature",values = myColors)
K <-  this_param_set$B.max

# LRR plot
g5A <- ggplot(temperature_data, aes(x = Days, y = interact_metric, colour = T.))+
  geom_line(size = 2)+
  colScale+
    ylim(-1.25,1.25)+
  geom_hline(yintercept = 0, lty = 2, size = 0.8)+
  xlab("Time (days)") + 
  theme(plot.title = element_text(hjust = 0.5))+
  ylab(expression(rho)) +
  theme_clean()+
  theme(axis.text.x = element_text(size = 14))+
  theme(axis.text.y = element_text(size = 14))+
  theme(axis.title.x = element_text(size = 18))+
  theme(axis.title.y = element_text(size = 18))


names(myColors) <- rev(levels(light_data$I.))
colScale <- scale_colour_manual(name = "Light",values = myColors)
K <-  this_param_set$B.max

g5B <- ggplot(light_data, aes(x = Days, y = interact_metric, colour = I.))+
  geom_line(size = 2)+
  colScale+
  geom_hline(yintercept = 0, lty = 2, size = 0.8)+
   ylim(-1.25,1.25)+
  xlab("Time (days)")  + 
  theme(plot.title = element_text(hjust = 0.5))+
  labs(color = "Light") +
  ylab(expression(rho)) +
  theme_clean()+
  theme(axis.text.x = element_text(size = 14))+
  theme(axis.text.y = element_text(size = 14))+
  theme(axis.title.x = element_text(size = 18))+
  theme(axis.title.y = element_text(size = 18))

})


# ---------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------
# Specify stressor scenarios to use across all models 
# ---------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------

# Generate the treatments we want to simulate
tempsens <- c(35, seq(36, 42, by = 2))
lightsens <- seq(200, 1000, by = 200)
templight <- expand.grid(T. = tempsens, I. = lightsens)

# -------------------------------------
# base model parameters 
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

this_param_set$T. <- templight$T.
this_param_set$I. <- templight$I.


# Base model plots --------------------------------------------------------

xout <- pmap(this_param_set, biomass.int.wrapper.cons.foodweb, model = "base") 

eval(make_plots_expr) # Evaluate the expression to generate plots

# bind plots to another name for when we recompute the expression
base_temp_plot <- g5A+ ggtitle("Base model")
base_light_plot <- g5B+ ggtitle("Base model")


# Holling type II model plots ---------------------------------------------

# Changes to parameter set
hol_params <- within(this_param_set,{
  h <- 1 # 
  })

# Run the model for each stressor scenario
xout <- pmap(hol_params, biomass.int.wrapper.cons.foodweb, model = "holling2") 

eval(make_plots_expr) # Evaluate the expression to generate plots

# bind plots to another name for when we recompute the expression
hol_temp_plot <- g5A+ ggtitle("Holling Type II model")
hol_light_plot <- g5B+ ggtitle("Holling Type II model")

# Gompertz model plots ----------------------------------------------------

# Changes to parameter set
gomp_params <- within(this_param_set,{
  
})

xout <- pmap(gomp_params, biomass.int.wrapper.cons.foodweb, model = "gompertz") 

eval(make_plots_expr) # Evaluate the expression to generate plots

# bind plots to another name for when we recompute the expression
gomp_temp_plot <- g5A+ ggtitle("Gompertz model")
gomp_light_plot <- g5B+ ggtitle("Gompertz model")

# Pmax dependent model plots ----------------------------------------------

# Changes to parameter set
pdep_params <- within(this_param_set,{
  
  psi <- tanh(1000/Ik)*Photo.Temp(T.= 35, PT.max, T.opt, T.max)/B.max

  })

xout <- pmap(pdep_params, biomass.int.wrapper.cons.foodweb, model = "pmaxdepend") 

eval(make_plots_expr) # Evaluate the expression to generate plots

# bind plots to another name for when we recompute the expression
pdep_temp_plot <- g5A+ ggtitle("Pmax model")
pdep_light_plot <- g5B + ggtitle("Pmax model")


# Summary of model plots --------------------------------------------------

# Compare temperature findings
(base_temp_plot + gomp_temp_plot) / ( pdep_temp_plot + hol_temp_plot)+plot_annotation(tag_levels = 'A')

ggsave(filename = "Plots/model_sens_consumer_temperature.png", width = 8, height = 6, units = c("in"), dpi = 300)

# Compare light findings
(base_light_plot + gomp_light_plot) / (pdep_light_plot + hol_light_plot)+plot_annotation(tag_levels = 'A')

ggsave(filename = "Plots/model_sens_consumer_light.png", width = 8, height = 6, units = c("in"), dpi = 300)


# Function to change the time limit to the first 50 days
ltrans <- function(x){
  x + xlim(c(0, 50))
}

## Transient dynamics of the interaction plot

# Compare temperature findings
(ltrans(base_temp_plot) + ltrans(gomp_temp_plot)) / 
  ( ltrans(pdep_temp_plot) + ltrans(hol_temp_plot))+plot_annotation(tag_levels = 'A')

ggsave(filename = "Plots/transient_model_sens_consumer_temperature.png", width = 8, height = 6, units = c("in"), dpi = 300)

# Compare light findings
(ltrans(base_light_plot) + ltrans(gomp_light_plot)) / 
  (ltrans(pdep_light_plot) + ltrans(hol_light_plot))+plot_annotation(tag_levels = 'A')

ggsave(filename = "Plots/transient_model_sens_consumer_light.png", width = 8, height = 6, units = c("in"), dpi = 300)
