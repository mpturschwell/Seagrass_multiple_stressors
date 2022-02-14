library(tidyverse)
library(RColorBrewer)
library(patchwork)
library(ggthemes)

source("Functions/conversion_functions.r")
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
f3A <- ggplot(experiment_data, aes(x = t, y = B, colour = Stressor))+
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

f3A
