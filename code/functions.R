library(markovchain)
library(msm)

force_of_infection = function(b, lambda, a, c, X, g, n){
  m=lambda/g
  (b * m * (a ^ 2) * c * X * exp(-g*n)) / (g + (a*c*X))
}

ac_calc = function(af, ap, rho, C, alpha, qt, qu, tau){
  a = af+ap
  D = C*rho + (C * (1-rho) * (qt)/(alpha*a + qt)) + ((1-C)*(qu / (a+qu)))
  delta = ((tau*C*rho) + (((1/qt) + tau)*(C * (1-rho) * (qt)/(alpha*a + qt))) +
   (((1/qu)+tau)*((1-C)*(qu / (a+qu))))) / D
  ac_inv = (delta * (D / (1-D))) + (1/a)
  return(1/ac_inv)
}


gc_calc = function(mu, qu, qt, gu, gt, gtau, tau, rho, C, pi.tau, pi.U, pi.T){
  return((pi.T * (mu*qt + gt)) + (pi.U * gu) + (pi.tau*gtau))
}

prob_untreated_house = function(mu, qu, qt, gu, gt, gtau, tau, rho, C){
  pi.tau = - (tau * qu * qt)/(qu * C * (rho-1) + qt*(C-tau*qu - 1))
  pi.U = (1-C) * pi.tau / (tau * qu)
  return(pi.U)
}

prob_treated_house = function( qu, qt, tau, rho, C){
  pi.tau = - (tau * qu * qt)/(qu * C * (rho-1) + qt*(C-tau*qu - 1))
  pi.U = (1-C) * pi.tau / (tau * qu)
  pi.t = (1-rho) * C * pi.tau / (tau * qt)
  return(pi.t)
}

prob_treated_house(qu=.1, qt=.1,tau=.1, rho=0, C=.3)
prob_treated_house(qu=.1, qt=.1, tau=.1, rho=.1, C=.3)
prob_treated_house( qu=.1, qt=.2, tau=.1, rho=0, C=.3)



parity = function(a, g){
  return(a / (a+g))
}

location_probability = function(af, ap, qu, tau, d, rho, qt, alpha, C) {
  tau.inv = 1/tau
  inf.matrix = matrix(c(
    -(d+qu*(1+d)), 0, d, qu, 0, qu*d, 0, 0, 0,
    
    ap, -(ap+qu*(ap+1)+2*d+qu*2*d), 2*d, qu*ap, qu, qu*2*d,0,0,0,
    
    af, ap, -(af+ap+qu*(af+ap+1)), qu*af, qu*ap, qu, 0,0,0,
    
    tau.inv*(1-C), 0, tau.inv*(1-C)*d, -(tau.inv*(1-C)*(1+d)+d+tau.inv*(1-rho)*C*(1+d)), 0, d, tau.inv*(1-rho)*C, 0, tau.inv*(1-rho)*C*d,
    
    tau.inv*(1-C)*ap, tau.inv*(1-C), tau.inv*(1-C)*2*d, ap, -(tau.inv*(1-C)*(ap+1)+ap+tau.inv*(1-rho)*C*(ap+1)+tau.inv*(1-C)*2*d+2*d+tau.inv*(1-rho)*C*2*d), 2*d, tau.inv*(1-rho)*C*ap, tau.inv*(1-rho)*C,tau.inv*(1-rho)*C*2*d,
    
    tau.inv*(1-C)*af, tau.inv*(1-C)*ap, tau.inv*(1-C), af, ap, -(tau.inv*(1-C)*(af+ap+1)+af+ap + tau.inv*(1-rho)*C*(af+1+ap)), tau.inv*(1-rho)*C*af, tau.inv*(1-rho)*C*ap, tau.inv*(1-rho)*C,
    
    0, 0, 0, qt, 0, qt*d, -(qt*(1+d)+d), 0, d,
    
    0,0,0, qt*ap*alpha, qt, qt*2*d, ap*alpha, -(qt*(ap*alpha+1)+ap*alpha+2*d+qt*d*2), 2*d,
    
    0,0,0,qt*af*alpha, qt*ap*alpha, qt, af*alpha,ap*alpha, -(qt*alpha*(af+ap)+qt+(af+ap)*alpha)
  ), nrow=9, ncol=9,byrow=T)
  prob_transition_matrix = MatrixExp(inf.matrix)
  markovchain_object = new("markovchain",
                           c("UF", "UP", "UN", "DF", "DP", "DN", "TF", "TP", "TN"),
              byrow=T, transitionMatrix = prob_transition_matrix)
  steadyStates= steadyStates(markovchain_object)
  if(nrow(steadyStates) == 1){
    return(steadyStates(markovchain_object))
  } else{
  return(steadyStates(markovchain_object)[2,])
  }
}

