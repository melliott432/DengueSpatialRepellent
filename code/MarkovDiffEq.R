library(deSolve)
library(markovchain)
library(msm)

feeding.markov.model = function(t, x, parms){
  
  UF = x[1]
  UP = x[2]
  UN = x[3]
  deltaF = x[4]
  deltaP = x[5]
  deltaN = x[6]
  TF = x[7]
  TP=x[8]
  TN=x[9]
  
  af = parms["af"]
  ap = parms["ap"]
  qu = parms["qu"]
  tau = parms["tau"]
  tau.inv = 1/tau
  d = parms["d"]
  rho=parms["rho"]
  qt = parms["qt"]
  alpha=parms["alpha"]
  C=params["C"]
  
  dUF = -(d+qu*(1+d))*UF + ap*UP + af*UN + (tau.inv*(1-C))*deltaF +
    tau.inv*(1-C)*ap*deltaP + tau.inv*(1-C)*af*deltaN
  dUP = -(ap+qu*(ap+1)+2*d+qu*2*d)*UP + ap*UN + tau.inv*(1-C)*deltaP +
    tau.inv*(1-C)*ap*deltaN
  
  dUN = d*UF +(d*2)*UP - (af+ap+qu*(af+ap+1))*UN+tau.inv*(1-C)*d*deltaF+
    tau.inv*(1-C)*(d*2)*deltaP + tau.inv*(1-C)*deltaN
  ddeltaF =  qu*UF + qu*ap*UP + qu*af*UN - (tau.inv*(1-C)*(1+d) + d +
                                              tau.inv*(1-rho)*C*(1+d))*deltaF+
    ap*deltaP+af*deltaN+qt*TF+qt*ap*alpha*TP + qt*af*alpha*TN
  
  ddeltaP = qu*UP + qu*ap*UN - (tau.inv*(1-C)*(ap+1) + ap +
    tau.inv*(1-rho)*C*(ap+1)+tau.inv*(1-C)*(d*2)+d*2+tau.inv*(1-rho)*C*(d*2))*deltaP +
    ap*deltaN+
    qt*TP + qt*ap*alpha*TN
  
  ddeltaN = qu*d*UF + qu*(d*2)*UP + qu*UN + d*deltaF +(d*2)*deltaP - 
    (tau.inv*(1-C)*(af+ap+1) + af +ap+tau.inv*(1-rho)*C*(af+ap+1))*deltaN +
    qt*d*TF + qt*(d*2)*TP + qt*TN
  
  dTF = tau.inv*(1-rho)*C*deltaF + tau.inv*(1-rho)*C*ap*deltaP + 
    tau.inv*(1-rho)*C*af*deltaN - (qt*(1+d) + d)*TF + ap*alpha*TP +
    af*alpha*TN
  
  dTP = tau.inv*(1-rho)*C*deltaP +
    tau.inv*(1-rho)*C*ap*deltaN-  
    (qt*(ap*alpha+1) + ap*alpha + (d*2) + qt*(d*2))*TP+
    ap*alpha*TN
  
  dTN = tau.inv * (1-rho) * C*d*deltaF +tau.inv*(1-rho)*C*(d*2)*deltaP +
    tau.inv*(1-rho)*C*deltaN+
    d*TF +(d*2)*TP - (qt*alpha*(af+ap)+qt +(af+ap)*alpha)*TN
  
  
  list(c(dUF, dUP, dUN, ddeltaF, ddeltaP, ddeltaN, dTF, dTP, dTN))
  
}

Time= seq(0,100, by=.1)
Init = c(.3, .3, .1, .1 ,.1, 0, .1,0,0)
qu=.3
tau=1
tau.inv = 1/tau
ap=0.3
af=.1
d=.3
rho = 0
qt=.2
alpha=.5
C=0.75

params=c(af=af, ap=ap, qu=qu, tau=tau, d=d, rho=rho, qt=qt, alpha=alpha, C=C)

out=as.data.frame(ode(Init, Time, feeding.markov.model, parms=params))

# sensitive to starting values
UF_prop = tail(out[,2],1)
UP_prop = tail(out[,3],1)
UN_prop = tail(out[,4],1)
deltaF_prop = tail(out[,5],1)
deltaP_prop = tail(out[,6],1)
delta_Nprop = tail(out[,7],1)
TF_prop= tail(out[,8],1)
TP_prop = tail(out[,9],1)
TN_prop = tail(out[,10],1)

UF_prop+UP_prop+UN_prop+deltaF_prop +deltaP_prop + delta_Nprop+
  TF_prop + TP_prop + TN_prop

plot(out[,2], type="l", ylim=c(0,1), main="Feeding Status- Indoor mosquitoes",
     xlab='Time', ylab="Proportion")
lines(out[,3], col="red")
lines(out[,4], col="blue")
legend("topright", c("Fully fed", "Partially fed", "Unfed"),
       col=c("black", "red", "blue"),
       lty=1)

# check function
c(UF_prop,UP_prop,UN_prop,deltaF_prop ,deltaP_prop , delta_Nprop,
  TF_prop , TP_prop , TN_prop)
location_probability(af, ap, qu, tau, d, rho, qt, alpha, C)
