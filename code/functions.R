library(markovchain)
library(msm)

force_of_infection = function(b, lambda, a, c, X, g, n){
  m=lambda/g
  (b * m * (a ^ 2) * c * X * exp(-g*n)) / (g + (a*c*X))
}

ac_calc = function(a, alpha, d, pi.UN, pi.TN, pi.tauN){
  au = a * d / (a + d)
  at = alpha * a * d / ((alpha*a) + d)
  
  return((pi.UN + pi.tauN)*au + pi.TN*at)
}


gc_calc = function(gu, g_mult, pi.tau, pi.U, pi.T){
  gt = gu * g_mult
  return(pi.T * gt + pi.U * gu + pi.tau*gu)
}

parity = function(a, g){
  return(a / (a+g))
}

location_probability = function(a, qu, tau, d, rho, q_mult, alpha, C) {
  inf.matrix = matrix(c(
    
    # untreated, blood fed
    -(d+qu), d, qu, 0, 0, 0,
    
    # untreated, not blood fed
    a, -(a + qu), 0, qu, 0, 0,
    
    # Transition, Blood Fed
    tau * (1-C), 0, -(d + tau*(1-C) + tau * C * (1-rho)), d, tau * C * (1-rho), 0,
    
    # Transition, not blood fed
    0, tau * (1-C), a, -(a + tau*(1-C) + tau*C*(1-rho)), 0, tau * C * (1-rho),
    
    # Treated, blood fed
    0, 0, q_mult*qu, 0, -(qu*q_mult + d), d,
    
    # Treated, not blood fed
    0, 0, 0, q_mult * qu, alpha * a, -(q_mult * qu + alpha*a)
    
  ), nrow=6, ncol=6,byrow=T)
  prob_transition_matrix = MatrixExp(inf.matrix)
  markovchain_object = new("markovchain",
                           c("UF", "UN", "DF",  "DN", "TF", "TN"),
                           byrow=T, transitionMatrix = prob_transition_matrix)
  steadyStates= steadyStates(markovchain_object)
  if(nrow(steadyStates) == 1){
    return(steadyStates(markovchain_object))
  } else{
    return(steadyStates(markovchain_object)[1,])
  }
}

