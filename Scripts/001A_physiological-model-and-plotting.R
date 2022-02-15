

library(purrr)
library(tidyverse)
library(ggthemes)
library(patchwork)

source("Functions/resp_function.R")
source("Functions/conversion_functions.R")
source("Functions/photo_light_shade_function.R")
source("Functions/photosyn_shade_function.R")
source("Functions/photo_temp_function.R")
source("Functions/net.growth_shade_function.R")

mytime <- format(Sys.time(), "%Y_%m_%d")

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


D_R <- function(T., I., B_range, T.control = 35, I.control = 1000, params = this_param_set){
  
  lwr_B <-  B_range[[1]]
  upr_B <-  B_range[[2]]
  df <- expand.grid(T. = c(T., T.control), I. = c(I., I.control), Ca = c(lwr_B, upr_B))
  
  params[["Ca"]] <- df$Ca
  params[["T."]] <- df$T.
  params[["I."]] <- df$I.
  
  net_growth <- unlist(pmap(params, Net.growth.shade))
  
  Y_x <- function(T_type, I_type, Ca_type){
    net_growth[ df$T. == ifelse(T_type == "C", T.control, T.) & 
                  df$I. == ifelse(I_type == "C", I.control, I.) & 
                  df$Ca == Ca_type]
  }
  
  D_R_lwr <- Y_x("C", "C", lwr_B) +  Y_x("S", "S", lwr_B) - 
    Y_x("S", "C", lwr_B) - Y_x("C", "S", lwr_B)
  D_R_upr <- Y_x("C", "C", upr_B) +  Y_x("S", "S", upr_B) - 
    Y_x("S", "C", upr_B) - Y_x("C", "S", upr_B)
  
  return(-c(D_R_lwr, D_R_upr))
  
}

B_range <- c(0, this_param_set[["B.max"]]*1.25) 

# Plot a line for ggplot
line_TI <- function(T., I.){
  geom_line(aes(x = B_range/this_param_set[["B.max"]]*100, y = D_R(T. = T., I. = I., B_range = B_range), 
                colour = paste0("Temp = ", T., ", Light = ", I.) ), size = 2)
}

colours <-  c("Temp = 42, Light = 200" = "#E31A1C",
              "Temp = 40, Light = 200" = "#FD8D3C",
              "Temp = 38, Light = 200" = "#FECC5C",
              "Temp = 36, Light = 200" = "#FFFFB2")

p1a <- ggplot() + 
  line_TI(42, 200) + 
  line_TI(40, 200) + 
  line_TI(38, 200) + 
  line_TI(36, 200) + 
  theme_bw() +
  geom_hline(yintercept = 0, lty = "dashed") + geom_vline(xintercept = c(0,100), lty = "dotted") +
  labs(y = expression(delta), color = "Treatment") +
  xlab(bquote('Percentage of carrying capacity')) +
  scale_color_manual(values = colours)+#,  guide = guide_legend(reverse = TRUE))+
  theme_clean()+
  theme(axis.text.x = element_text(size = 14))+
  theme(axis.text.y = element_text(size = 14))+
  theme(axis.title.x = element_text(size = 14))+
  theme(axis.title.y = element_text(size = 18))+
  ylim(c(-0.08,0.02))

p1a

colours2 <-  c("Temp = 42, Light = 200" = "#E31A1C",
               "Temp = 42, Light = 400" =  "#FD8D3C",
               "Temp = 42, Light = 600" =  "#FECC5C",
               "Temp = 42, Light = 800" =  "#FFFFB2")

p1b <- ggplot() + 
  line_TI(42, 200) + 
  line_TI(42, 400) + 
  line_TI(42, 600) + 
  line_TI(42, 800) + 
  theme_bw() +
  geom_hline(yintercept = 0, lty = "dashed") + geom_vline(xintercept = c(0,100), lty = "dotted") +
  labs(y = expression(delta), color = "Treatment") +
  xlab(bquote('Percentage of carrying capacity')) +
  scale_color_manual(values = colours2)+#,  guide = guide_legend(reverse = TRUE))+
  theme_clean()+
  theme(axis.text.x = element_text(size = 14))+
  theme(axis.text.y = element_text(size = 14))+
  theme(axis.title.x = element_text(size = 18))+
  theme(axis.title.y = element_text(size = 18))

p1b

p1 <- p1a / p1b
p1 <- p1 + plot_annotation(tag_levels = 'A')
p1

ggsave(path = "Plots", filename = paste0(mytime, "_Figure2.tiff"), p1, width = 8, height = 8, units = c("in"), dpi = 300)


