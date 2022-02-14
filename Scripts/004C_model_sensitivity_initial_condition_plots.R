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
source("Functions/conversion_functions.r")
source("Functions/get_stressor_interaction.R")
source("Functions/photo_temp_function.R")
source("Functions/resp_function.R")
source("Functions/photo_light_shade_function.R")
source("Functions/photosyn_shade_function.R")
source("Functions/biomass.integration_shade_function_consumer_foodweb.R")
source("Functions/dBdt_shade_function_consumer_foodweb.R")
source("Functions/attack_mort_temp_scaling_function.R")
source("Functions/alternative_consumer_models_dBdt.R")


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
tempsens <- c(35, 42)
lightsens <- c(1000, 200)
Bsens <- c(145.75, 435.25)
Xsens <- c(33.25, 97.75)
templight <- expand.grid(T. = tempsens, I. = lightsens, B_init = Bsens, X_init = Xsens)

# Remove the treatments we don't need
#templight <- templight[templight$T. %in% c(35, 42.58) | templight$I. %in% c(254.55, 1000),] 

this_param_set$T. <- templight$T.
this_param_set$I. <- templight$I.
this_param_set$B_init <- templight$B_init
this_param_set$X_init <- templight$X_init


get_model_IC_data <- function(param_set, model){
  
  xout <- pmap(param_set, biomass.int.wrapper.cons.foodweb, model = model) 
  
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
    Days = (time-1)/24
    B_init <- factor(init)
    })
  
  return(out)
}


# Get the IC changing dataframe from each model and bind them
all_mod_IC_df <- c("base", "gompertz", "pmaxdepend", "holling2")

base_df <- this_param_set %>% get_model_IC_data(model = "base") %>%
  mutate(model = "base")
gompertz_df <- this_param_set %>% get_model_IC_data(model = "gompertz") %>% 
  mutate(model = "gompertz")
pmaxdepend_df <- this_param_set %>% 
  within({psi = tanh(1000/Ik)*PT.max*(T.max-35)/(T.max - T.opt) *
           (35/T.opt)^(T.opt/(T.max - T.opt))/B.max}) %>% 
  get_model_IC_data(model = "pmaxdepend") %>% mutate(model = "pmaxdepend")
holling2_df <- this_param_set %>% 
  within({
    v <- 0.0275
    h <- 3
    }) %>% 
  get_model_IC_data(model = "holling2") %>% mutate(model = "holling2")

all_mod_IC_df <- bind_rows(base_df, gompertz_df, pmaxdepend_df, holling2_df)

K <-  this_param_set$B.max


library(RColorBrewer)
myColors <- brewer.pal(4,"Dark2")
names(myColors) <- levels(all_mod_IC_df$init)
colScale <- scale_colour_manual(name = "Initial Conditions",values = myColors)

# LRR plot
g1 <- ggplot(all_mod_IC_df, aes(x = Days, y = interact_metric, colour = init))+
  geom_line(size = 2)+
  colScale+
  #ylim(-1.85,0.25)+
  geom_hline(yintercept = 0, lty = 2, size = 0.8)+
  xlab("Time (days)") + 
  theme(plot.title = element_text(hjust = 0.5))+
  ylab(expression(rho)) +
  theme_clean()+
  theme(axis.text.x = element_text(size = 14))+
  theme(axis.text.y = element_text(size = 14))+
  theme(axis.title.x = element_text(size = 18))+
  theme(axis.title.y = element_text(size = 18))+
  facet_wrap(.~ model) + xlim(c(0, 250))
g1


ggsave(path = "Plots", filename = "model_sens_initial_conditions_1.png", g1, width = 30, height = 25, units = c("cm"), dpi = 300)
