dBdt.gompertz <- function(t, B, X, params){
  # Model with gompertz growth shading
  with(as.list(params), { 
    P <- photosyn.shade(T., I., PT.max, PL.max, T.opt, T.max, Ik, Ca = B, B.max)
    P <- Photo.Temp(T., PT.max, T.opt, T.max)*tanh(I./Ik)*log(B.max/B) 
    R <- Resp.temp(R.max, T., RT.opt, RT.max)
    S <- Attack.Mortal.Temp.Scaling(T.)
    dB <- (P - R - M - (a*S*X))*B
    dX <- (a*S*c*X*B)-(v*S*X)
    c(dB,dX)
  })
}

dBdt.pmaxdepend <- function(t, B, X, params){
  # Model with zero photosynthesis dependent on Pmax
  with(as.list(params), { 
    P <- Photo.Temp(T., PT.max, T.opt, T.max)*tanh(I./Ik) - psi*B 
    R <- Resp.temp(R.max, T., RT.opt, RT.max)
    S <- Attack.Mortal.Temp.Scaling(T.)
    dB <- (P - R - M - (a*S*X))*B
    dX <- (a*S*c*X*B)-(v*S*X)
    c(dB,dX)
  })
}

dBdt.holling2 <- function(t, B, X, params){
  # Model with a holling type II functional response for predator
  with(as.list(params), { 
    P <- photosyn.shade(T., I., PT.max, PL.max, T.opt, T.max, Ik, Ca = B, B.max)
    R <- Resp.temp(R.max, T., RT.opt, RT.max)
    S <- Attack.Mortal.Temp.Scaling(T.)
    dB <- (P - R - M)*B - (a*S*B)/(1 + a*S*h*B)*X
    dX <- (a*S*B)/(1 + a*S*h*B)*c*X-(v*S*X)
    c(dB,dX)
  })
}
