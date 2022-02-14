
# Max's stressor interaction function

get_stressor_interaction <- function(data, temp, light, T_control = 35, I_control = 1000){
  
  if ("biomass" %in% names(data)) {
    names(data)[names(data) == "biomass"] <- "B"
  }
  
  control <- dplyr::filter(data, T. == T_control, I. == I_control)
  temp_stress <- dplyr::filter(data, T. == T_control, I. == light)
  light_stress <- dplyr::filter(data, T. == temp, I. == I_control)
  both_stress <- dplyr::filter(data, T. == temp, I. == light)
  
  
  if ("netgrowth" %in% names(data)){
    # in this case we cannot using the multiplicative metric doesn't make sense so we use the additive metric
    both_stress$interact_metric <-  (both_stress$netgrowth - control$netgrowth) - 
      ((light_stress$netgrowth) + (temp_stress$netgrowth) - 2*control$netgrowth)
    
  } else {
    # We use the multiplicative metric I_R
    both_stress$interact_metric <-  (log(both_stress$B) - log(control$B)) - 
      (log(light_stress$B) + log(temp_stress$B) - 2*log(control$B))
    
  }
  
  # Flip the sign of the interaction metric - at the request of the reviewer
  both_stress$interact_metric <- -both_stress$interact_metric
  
  return(both_stress)
  
}