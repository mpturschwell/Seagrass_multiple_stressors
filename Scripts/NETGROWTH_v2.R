#----------------------------------------------#
# Sensitivty analysis for STRESSNET MODEL 

#----------------------------------------------#

# setwd and load libraries 
#setwd("/Users/s2904436/Dropbox/Multi-stressor/Code") # I just use the R-Project

library(purrr)
library(ggplot2)
library(dplyr)
library(patchwork)


# source functions 
source("Functions/resp_function.R")
source("Functions/MC_conversion_functions.R")
source("Functions/photo_light_shade_function.R")
source("Functions/photosyn_shade_function.R")
source("Functions/photo_temp_function.R")
source("Functions/net.growth_shade_function.R")
source("Functions/get_stressor_interaction.R")

# -------------------------------------
# Input required parameters 
# -------------------------------------

this_param_set <- list(
  Ca = 1,
  B.max = 667, 
  T. = 30,      # Temp                          degrees C
  I. = 319,     # Irradiance                    mmol m-2 second-1
  PT.max = mg.C.per.h_to_g.dry.wt.per.d(4.7),   # Max gross production (temp)   mg C g-1 dry wt h-1
  PL.max = umol_to_mg(422), # Max gross production (light)  mmol O2 g dry wt-1 hr-1
  T.opt = 34.9, # optimum temperature           degrees C
  T.max = 44.5, # Maximum temperature           degrees C
  Ik = 319,     # Saturation irradiance         mmol m-2 second-1
  R.max = mg.C.per.h_to_g.dry.wt.per.d(1.1),  # Maximum respiration           mg C g-1 dry wt h-1
  RT.opt = 39.1,  # Resp optimum temperature      degrees C
  RT.max = 45.6 # Resp Maximum temperature      degrees C
)
# -------------------------------------
# create parameter sequence to run model 
# -------------------------------------

Ca <- seq(30, 667, length.out = 100)                      # 10, 50 and 90% carrying capacity

# Generate the treatments we want to simulate
tempsens <- c(35, seq(36, 42, by = 2))
lightsens <- seq(200, 1000, by = 200)
templight <- expand.grid(T. = tempsens, I. = lightsens, Ca = Ca)

this_param_set$T. <- templight$T.
this_param_set$I. <- templight$I.
this_param_set$Ca <- templight$Ca

# --------------------------------------------------
# Run pmap and net growth function on all combinations 
# --------------------------------------------------
xout <- pmap(this_param_set, Net.growth.shade) %>% 
  unlist()

templight$netgrowth <- xout


# ----------------------------------------
# Plot outputs 
# --------------------------------------
templight$T. <- as.factor(templight$T.)
templight$I. <- as.factor(templight$I.)

# make labels 
control <- templight %>% 
  dplyr::filter(T. == "35", I. == "1000") %>% 
  mutate(Stressor = "Control")

Light <- templight %>% 
  dplyr::filter(T. == "35", I. == "200") %>% 
  mutate(Stressor = "Light Stress")

Temp <- templight %>% 
  dplyr::filter(T. == "42", I. == "1000") %>% 
  mutate(Stressor = "Temp Stress")

stress <- templight %>% 
  dplyr::filter(T. == "42", I. == "200") %>% 
  mutate(Stressor = "Temp + Light Stress")

Out <- rbind(control, Light, Temp, stress)
Out$Stressor <- as.factor(Out$Stressor)

Out$Stressor <- ordered(Out$Stressor, levels = c("Control", "Light Stress", "Temp Stress", "Temp + Light Stress"))

rescale <- function(x) (x-min(x))/(max(x) - min(x)) * 100
Out$Ca <- rescale(Out$Ca)


library(RColorBrewer)
myColors <- brewer.pal(4,"YlOrRd")
names(myColors) <- levels(Out$Stressor)
colScale <- scale_colour_manual(name = "Stressor",values = myColors)

g1A <- ggplot(Out, aes(x = Ca, y = netgrowth, colour = Stressor))+
  geom_line(size = 2)+
  #colScale+
  xlab("Carrying capacity (%)") + labs(colour = "Stressor") + 
  theme(plot.title = element_text(hjust = 0.5))+
  ylab(bquote('Net production rate ('*B^-1~d^-1*')')) +
  theme_bw()
g1A




###################
# -----------------------------------------------------------------------------------
#      Figure 3 -  Magnitude variation 
# -----------------------------------------------------------------------------------

light_for_temp <- 200
temp_for_light <- 42

templight
temperature_data <- rbind(get_stressor_interaction(data = templight, temp = 36, light = light_for_temp),
                          get_stressor_interaction(data = templight, temp = 38, light = light_for_temp),
                          get_stressor_interaction(data = templight, temp = 40, light = light_for_temp),
                          get_stressor_interaction(data = templight, temp = 42, light = light_for_temp)) %>% within({
                            T. <- factor(T., levels = c(42, 40, 38, 36))
                          })

light_data <- rbind(get_stressor_interaction(data = templight, temp = temp_for_light, light = 200),
                    get_stressor_interaction(data = templight, temp = temp_for_light, light = 400),
                    get_stressor_interaction(data = templight, temp = temp_for_light, light = 600),
                    get_stressor_interaction(data = templight, temp = temp_for_light, light = 800)) %>% within({
                      I. <- factor(I., levels = c(200, 400, 600, 800))
                    })


library(RColorBrewer)

names(myColors) <- rev(levels(temperature_data$T.))
colScale <- scale_colour_manual(name = "Temperature",values = myColors)

# LRR plot
g3A <- ggplot(temperature_data, aes(x = Ca, y = interact_metric, colour = T.))+
  geom_line(size = 2)+
  colScale+
  #  ylim(-1.25,1.25)+
  geom_hline(yintercept = 0, lty = 2, size = 0.8)+
  xlab("Carrying Capacity %") + 
  theme(plot.title = element_text(hjust = 0.5))+
  ylab(expression(D[R])) +
  theme_bw()
g3A


names(myColors) <- rev(levels(light_data$I.))
colScale <- scale_colour_manual(name = "Light",values = myColors)

g3B <- ggplot(light_data, aes(x = Ca, y = interact_metric, colour = I.))+
  geom_line(size = 2)+
  colScale+
  #  ylim(-1.25,1.25)+
  geom_hline(yintercept = 0, lty = 2, size = 0.8)+
  xlab("Carrying Capacity %") + 
  theme(plot.title = element_text(hjust = 0.5))+
  ylab(expression(D[R])) +
  theme_bw()
g3B

fig5 <- g3A/ g3B
fig5 <-fig5 + plot_annotation(tag_levels = 'A')
fig5
ggsave(path = "Plots", filename = "2021-06-15_NET_productivity_var_magnitude_TEMP-LIGHT.png", fig5, width = 10, height = 10, units = c("in"), dpi = 300)

