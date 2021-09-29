library(deSolve)


# define the gradient function
location_model = function(t, y, parms){
  # Pull state variables from y vector
  Delta=y[1]
  U=y[2]
  T=y[3]
  
  # Pull parameter values
  tau = parms["tau"]
  rho= parms["rho"]
  qu= parms["qu"]
  qt= parms["qt"]
  
  tau.inv = 1/tau
  
  # equations
  dDelta = tau.inv * (rho*C - 1)*Delta + qu*U + qt*T
  dU = -qu*U+tau.inv*(1-C)*Delta
  dT = -qt*T+tau.inv*(1-rho)*C*Delta
  
  # return list of gradients
  return(list(c(dDelta, dU, dT)))
}

times=seq(0, 10000, by=1/10)
parms = unlist(prior.samps[1,])

start=c(Delta=.3, U=.3, T=.4)

# integrate
out=as.data.frame(ode(y=start, t=times, parms=parms, func=location_model))

prop_transit= tail(out$Delta, 1)
prop_untreated = tail(out$U, 1)
prop_treated = tail(out$T, 1)

tau = parms["tau"]
rho= parms["rho"]
qu= parms["qu"]
qt= parms["qt"]

pi.U = (1-C) * pi.tau / (tau * qu)

# verified

