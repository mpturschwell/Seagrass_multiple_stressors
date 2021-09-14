# ------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------
#
#   Sensitivity for BIOMASS FOODWEB MODEL 
#
# ------------------------------------------------------------------------------------------------------
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
myColors <- brewer.pal(4,"YlOrRd")
names(myColors) <- levels(Out$Stressor)
colScale <- scale_colour_manual(name = "Stressor",values = myColors)

K <-  this_param_set$B.max
g1C <- ggplot(Out, aes(x = time*this_param_set[["dt"]], y = biomass, colour = Stressor))+
  geom_line(size = 2)+
  #colScale+
  xlab("Time (days)") + labs(colour = "Stressor") + 
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

# YOU MUST RUN Figure1B in the BIOMASS_v2.R script first
Figure2 <- g1B/g1C
Figure2 <- Figure2 + plot_annotation(tag_levels = 'A')
Figure2
ggsave(path = "Plots", filename = paste0(mytime, "_Figure3_250d.png"), Figure2 , width = 8, height = 8, units = c("in"), dpi = 300)




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


names(myColors) <- rev(levels(temperature_data$T.))
colScale <- scale_colour_manual(name = "Temperature",values = myColors)
K <-  this_param_set$B.max

# LRR plot
g5A <- ggplot(temperature_data, aes(x = Days, y = interact_metric, colour = T.))+
  geom_line(size = 2)+
  colScale+
  #  ylim(-1.25,1.25)+
  geom_hline(yintercept = 0, lty = 2, size = 0.8)+
  xlab("Time (days)") + 
  theme(plot.title = element_text(hjust = 0.5))+
  ylab(expression(I[R])) +
  theme_clean()+
  theme(axis.text.x = element_text(size = 14))+
  theme(axis.text.y = element_text(size = 14))+
  theme(axis.title.x = element_text(size = 14))+
  theme(axis.title.y = element_text(size = 14))
g5A


names(myColors) <- rev(levels(light_data$I.))
colScale <- scale_colour_manual(name = "Light",values = myColors)
K <-  this_param_set$B.max

g5B <- ggplot(light_data, aes(x = Days, y = interact_metric, colour = I.))+
  geom_line(size = 2)+
  colScale+
  geom_hline(yintercept = 0, lty = 2, size = 0.8)+
  # ylim(-1.25,1.25)+
  xlab("Time (days)")  + 
  theme(plot.title = element_text(hjust = 0.5))+
  ylab(expression(I[R])) + labs(color = "Light") +
  theme_clean()+
  theme(axis.text.x = element_text(size = 14))+
  theme(axis.text.y = element_text(size = 14))+
  theme(axis.title.x = element_text(size = 14))+
  theme(axis.title.y = element_text(size = 14))
g5B

fig5 <- g5A/ g5B
fig5 <-fig5 + plot_annotation(tag_levels = 'A')
fig5
ggsave(path = "Plots", filename = paste0(mytime, "_Figure5_FOODWEB_TEMP-LIGHT.png"), fig5, width = 10, height = 10, units = c("in"), dpi = 300)

