
B_t <- function(params, t){
  # Compute the analytical solution at each time t based on the parameters
  # As per "Population_model_solution.docx" by Max Campbell
  
  with(params, {
  
    phi <-  tanh(I./Ik)*PT.max*(T.max-T.)/(T.max - T.opt) * (T./T.opt)^(T.opt/(T.max - T.opt))
    R <- R.max*(RT.max-T.)/(RT.max-RT.opt) * (T./RT.opt) ^ (RT.opt/(RT.max-RT.opt)) 
    g <- phi - R - M
    alph <- -phi/B.max
    F. <-  B_init/(alph*B_init + g)
  
    return(g*F.*exp(g*t)/(1 - alph*F.*exp(g*t)))
  })
}


# Other plausible population models ---------------------------------------

 
B_t_gompertz <- function(params, t){
  # Compute at time t using the analytical solution for the gompertz model
  # see the doc ""
  
  with(params, {
    
    P.max <-  tanh(I./Ik)*PT.max*(T.max-T.)/(T.max - T.opt) * (T./T.opt)^(T.opt/(T.max - T.opt))
    R <- R.max*(RT.max-T.)/(RT.max-RT.opt) * (T./RT.opt) ^ (RT.opt/(RT.max-RT.opt)) 
    output <- B.max*exp(-(R + M)/P.max + (log(B_init/B.max) + (R + M)/P.max)*exp(-P.max*t))
    
    return(output)
  })
}

B_t_pmaxdepend <- function(params, t){
  # Compute at time t using the analytical solution for the gompertz model
  # see the doc ""
  
  with(params, {
    
    P.max <-  tanh(I./Ik)*PT.max*(T.max-T.)/(T.max - T.opt) * (T./T.opt)^(T.opt/(T.max - T.opt))
    R <- R.max*(RT.max-T.)/(RT.max-RT.opt) * (T./RT.opt) ^ (RT.opt/(RT.max-RT.opt)) 
    
    output <- ((P.max - R - M)/psi)/(1-(1-(P.max - R - M)/(psi*B_init))*exp(-(P.max - R - M)*t))
      
    return(output)
    
  })
}

