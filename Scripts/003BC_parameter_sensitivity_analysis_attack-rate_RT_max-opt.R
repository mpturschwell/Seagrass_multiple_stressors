# ------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------
#
#   Sensitivity for BIOMASS and Consumer MODELs 
#
# ------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------


library(tidyverse)
library(patchwork)
library(gridExtra)
library(RColorBrewer)


# source functions 
source("Functions/photo_temp_function.R")
source("Functions/resp_function.R")
source("Functions/photo_light_shade_function.R")
source("Functions/photosyn_shade_function.R")
source("Functions/biomass.integration_shade_function_consumer_foodweb.R")
source("Functions/dBdt_shade_function_consumer_foodweb.R")
source("Functions/attack_mort_temp_scaling_function.R")
source("Functions/population_model_analytical.R")
source("Functions/conversion_functions.R")
source("Functions/get_stressor_interaction.R")


# -------------------------------------
# set base parameters 
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


# Choose the sensitivity analyses to perform ------------------------------

# How many temps between T.opt and T.max
n_temps <- 5

# Choose here
sens_params <- list(
T.opt = c(20, 25, 30, 40),
T.max = c(35, 40, 50),
RT.opt = c(25, 30, 35, 45),
RT.max = c(40, 45, 50),
a = c(0.00001, 0.0001, 0.01, 0.1)
)

sens_params <- unlist(sens_params)

IR_plot_list <-  list(biomass_model = NULL, consumer_model = NULL)


# For each of the sensitivity scenarios
for (i in seq_along(sens_params)){

# Process the different sensitivity scenarios -----------------------------


param_changed <- substr(names(sens_params)[[i]], start = 1, stop = nchar(names(sens_params)[[i]]) - 1)
sens_scenario <- within(this_param_set, assign(x = param_changed, value = sens_params[[i]]))

# Generate the treatments we want to simulate
min_val <- ifelse(param_changed %in% c("RT.max", "RT.opt"), 
                  sens_scenario$RT.opt, sens_scenario$T.opt)
max_val <- min(sens_scenario$T.max, sens_scenario$RT.max) # to ensure T <= T.max and T <= RT.max
tempsens <- unique(c(seq(min_val, max_val, length.out = min(n_temps, max(floor(max_val-min_val), 1))), 
                     sens_scenario$T.opt)) 
tempsens_ext <- unique(tempsens, sens_scenario$T.opt) # Add this so we can compute the I_R
lightsens <- c(200, 1000)
templight <- expand.grid(T. = tempsens_ext, I. = lightsens)

sens_scenario$T. <- templight$T.
sens_scenario$I. <- templight$I.


# Run the consumer model for the parameters ---------------------------------------

xout <- pmap(sens_scenario, biomass.int.wrapper.cons.foodweb) 

# extract seagrass biomass vectors only 
xout2 <- lapply(xout, function(x) x[1][[1]])

xout3 <- xout2 %>%
  do.call("rbind", .) %>%
  data.frame() %>%
  cbind(templight) %>%
  tidyr::gather(time, biomass, -T., -I.) %>%
  dplyr::mutate(time = as.numeric(substring(time, 2)))

light_for_temp <- 200

consumer_data <- lapply(tempsens, function (x){ 
  get_stressor_interaction(data = xout3, temp = x, light = light_for_temp, 
                           T_control = sens_scenario[["T.opt"]])
    }) %>% do.call("rbind", .) %>% within({
      Days <-  (time-1)*this_param_set[["dt"]]
      T. <- factor(T., levels = rev(tempsens), labels = round(rev(tempsens),1))
    })



# Run the biomass model for the parameters --------------------------------

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

experiment_data2 <- compare_param(params = sens_scenario, t = seq(0, 365, by = 1),
                                  T. = templight$T., I. =templight$I.) %>% 
  do.call("rbind", .)


biomass_data <- lapply(tempsens, function (x){ 
  get_stressor_interaction(data = experiment_data2, temp = x, light = light_for_temp, 
                           T_control = sens_scenario[["T.opt"]])
}) %>% do.call("rbind", .) %>% within({
  Days <-  t
  T. <- factor(T., levels = rev(tempsens), labels = round(rev(tempsens),1))
})


# Interaction plots --------------------------------------------------------

temp_dif <- round(tempsens - sens_scenario[["T.opt"]]) + 26
cool <-  rainbow(25, start=rgb2hsv(col2rgb('cyan'))[1], end=rgb2hsv(col2rgb('blue'))[1])
warm <-  rainbow(25, start=rgb2hsv(col2rgb('red'))[1], end=rgb2hsv(col2rgb('yellow'))[1])
cols <- c(rev(cool), "grey", rev(warm))

myColors <- cols[temp_dif]
names(myColors) <- round(tempsens,1)
colScale <- scale_colour_manual(name = "Temperature",values = myColors)


# IR plot

g1 <- ggplot(biomass_data, aes(x = Days, y = interact_metric, colour = T.))+
  geom_line(size = 2)+
  colScale+
  geom_hline(yintercept = 0, lty = 2, size = 0.8)+
  # ylim(-1.25,1.25)+
  xlab("Time (days)")  + 
  theme(plot.title = element_text(hjust = 0.5))+
  ylab(expression(I[R])) + labs(color = "Light") +
  theme_bw() + ylim(-1.5, 1.5)

g2 <- ggplot(consumer_data, aes(x = Days, y = interact_metric, 
                                 colour = T.))+
  geom_line(size = 2)+
  colScale+
  geom_hline(yintercept = 0, lty = 2, size = 0.8)+
  xlab("Time (days)") + 
  theme(plot.title = element_text(hjust = 0.5))+
  ylab(expression(rho)) +
  theme_bw() + ylim(-1.5, 1.5)


IR_plot_list[["biomass_model"]][[i]] <- g1
IR_plot_list[["consumer_model"]][[i]] <- g2


}

bp <- function(x){ 
  IR_plot_list[["biomass_model"]][[x]] + 
    ggtitle(
      paste0(substr(names(sens_params)[[x]], 
                   start = 1, stop = nchar(names(sens_params)[[x]]) - 1), 
      " = ", sens_params[[x]])
    )
}
cp <- function(x) { 
  IR_plot_list[["consumer_model"]][[x]] + 
    ggtitle(
      paste0(substr(names(sens_params)[[x]], 
                    start = 1, stop = nchar(names(sens_params)[[x]]) - 1), 
             " = ", sens_params[[x]])
    )
}

fig1 <- (bp(1) + bp(2) + bp(3))/ (bp(4) + bp(5) + bp(6)) /
  (bp(7) + bp(8) + bp(9)) / (bp(10) + bp(11) + bp(12)) / 
  (bp(13) + bp(14) + plot_spacer())
#fig1 <-fig1 + plot_annotation(tag_levels = 'A')
fig1
ggsave(path = "Plots", filename = "2021-06-15_sens_temps_biomass.png", fig1, 
       width = 23, height = 30, units = c("cm"), dpi = 300)

fig2 <- (cp(1) + cp(2) + cp(3))/ (cp(4) + cp(5) + cp(6)) /
  (cp(7) + cp(8) + cp(9)) / (cp(10) + cp(11) + cp(12)) / 
  (cp(13) + cp(14) + plot_spacer())
#fig2 <-fig2 + plot_annotation(tag_levels = 'A')
fig2
ggsave(path = "Plots", filename = "2021-06-15_sens_temps_consumer.png", fig2, 
       width = 23, height = 30, units = c("cm"), dpi = 300)

fig3 <- (cp(15) + ylim(c(-1.5, 1.5)) + cp(16) + ylim(c(-1.5, 1.5)))/(cp(17) + ylim(c(-4, 4)) + cp(18) + ylim(c(-32, 20)))
#fig3 <-fig3 + plot_annotation(tag_levels = 'A')
fig3
ggsave(path = "Plots", filename = "2021-06-15_sens_attack_consumer.png", fig3, 
       width = 20, height = 15, units = c("cm"), dpi = 300)

# ggsave(path = "Plots", filename = "2021-06-15_sensi.png", fig5, width = 10, height = 10, units = c("in"), dpi = 300)

