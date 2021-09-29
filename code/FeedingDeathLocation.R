library(deSolve)
library(tidyverse)

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
  D = x[10]
  
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
  gu = params["gu"]
  gt = params["gt"]
  gtau = params["gtau"]
  
  dUF = -(d+qu*(1+d)+gu)*UF + ap*UP + af*UN + (tau.inv*(1-C))*deltaF +
    tau.inv*(1-C)*ap*deltaP + tau.inv*(1-C)*af*deltaN
  dUP = -(ap+qu*(ap+1)+gu)*UP + ap*UN + tau.inv*(1-C)*deltaP +
    tau.inv*(1-C)*ap*deltaN
  
  dUN = d*UF - (af+ap+qu*(af+ap+1)+gu)*UN+tau.inv*(1-C)*d*deltaF+
    tau.inv*(1-C)*deltaN
  ddeltaF =  qu*UF + qu*ap*UP + qu*af*UN - (tau.inv*(1-C)*(1+d) + d +
                                              tau.inv*(1-rho)*C*(1+d)+gtau)*deltaF+
    ap*deltaP+af*deltaN+qt*TF+qt*ap*alpha*TP + qt*af*alpha*TN
  
  ddeltaP = qu*UP + qu*ap*UN - (tau.inv*(1-C)*(ap+1) + ap +
                                  tau.inv*(1-rho)*C*(ap+1)+gtau)*deltaP + ap*deltaN+
    qt*TP + qt*ap*alpha*TN
  
  ddeltaN = qu*d*UF + qu*UN + d*deltaF - 
    (tau.inv*(1-C)*(af+ap+1) + af +ap+tau.inv*(1-rho)*C*(af+ap+1)+gtau)*deltaN +
    qt*d*TF + qt*TN
  
  dTF = tau.inv*(1-rho)*C*deltaF + tau.inv*(1-rho)*C*ap*deltaP + 
    tau.inv*(1-rho)*C*af*deltaN - (qt*(1+d) + d + gt)*TF + ap*alpha*TP +
    af*alpha*TN
  
  dTP = tau.inv*(1-rho)*C*deltaP +
    tau.inv*(1-rho)*C*ap*deltaN-  
    (qt*(ap*alpha+1) + ap*alpha + gt)*TP+
    ap*alpha*TN
  
  dTN = tau.inv * (1-rho) * C*d*deltaF + tau.inv*(1-rho)*C*deltaN+
    d*TF - (qt*alpha*(af+ap)+qt +(af+ap)*alpha + gt)*TN
  
  dD = gu*(UF+UP+UN) + gtau * (deltaF+deltaP+deltaN) + gt * (TF+TP+TN)
  
  
  list(c(dUF, dUP, dUN, ddeltaF, ddeltaP, ddeltaN, dTF, dTP, dTN, dD))
  
}

# test it out

Time= seq(0,40, by=.1)
Init = c(0,0,0,0,0,1,0,0,0,0)
qu=.3
tau=1
tau.inv = 1/tau
ap=0.2
af=.1
d=.1
rho = .2
qt=.2
alpha=.9
C=0
gu = .1
gtau=.1
gt=.2

params=c(af=af, ap=ap, qu=qu, tau=tau, d=d, rho=rho, qt=qt,
         alpha=alpha, C=C, gu=gu, gtau=gtau, gt=gt)

out=as.data.frame(ode(Init, Time, feeding.markov.model, parms=params))

#### All starting in transition unfed ####
plot(out[,2]~ out[,1], type="l", col="blue",lty="dashed",
     ylim=c(0,1), main="Compartment-specific death, 0 coverage",
     xlab="Time", ylab="Proportion Mosquitoes")
lines(out[,3]~ out[,1], type="l", col="red", lty="dashed")
lines(out[,4]~ out[,1], type="l", col="green", lty="dashed")
lines(out[,5]~ out[,1], type="l", col="blue")
lines(out[,6]~ out[,1], type="l", col="red")
lines(out[,7]~ out[,1], type="l", col="green")
lines(out[,11]~ out[,1], type="l", col="purple")

legend("right", c("Untreated House", "Transition", "Unfed", "Part-fed",
                  "Fully fed", "Dead"),
       lty=c(2, 1, 1, 1, 1, 1),
       col=c("black", "black", "green", "red", "blue", "purple"))

#### with treatment ####
qu=.3
tau=1
tau.inv = 1/tau
ap=0.2
af=.1
d=.1
rho = .2
qt=.2
alpha=.9
C=.4
gu = .1
gtau=.1
gt=.2

params=c(af=af, ap=ap, qu=qu, tau=tau, d=d, rho=rho, qt=qt,
         alpha=alpha, C=C, gu=gu, gtau=gtau, gt=gt)

out=as.data.frame(ode(Init, Time, feeding.markov.model, parms=params))

plot(out[,2]~ out[,1], type="l", col="blue",lty="dashed",
     ylim=c(0,1), main="Compartment-specific death, 40% coverage",
     xlab="Time", ylab="Proportion Mosquitoes")
lines(out[,3]~ out[,1], type="l", col="red", lty="dashed")
lines(out[,4]~ out[,1], type="l", col="green", lty="dashed")
lines(out[,8]~ out[,1], type="l", col="blue")
lines(out[,9]~ out[,1], type="l", col="red")
lines(out[,10]~ out[,1], type="l", col="green")
lines(out[,11]~ out[,1], type="l", col="purple")

legend("right", c("Untreated House", "Treated House", "Unfed", "Part-fed",
                  "Fully fed", "Dead"),
       lty=c(2, 1, 1, 1, 1, 1),
       col=c("black", "black", "green", "red", "blue", "purple"))


#### with an emergence rate ####


  

params=c(af=.2, ap=0, qu=.1, tau=.1, d=.3, rho=.3, qt=.2,
         alpha=.9, C=.5, gu=.1, gtau=.3, gt=.2)
out=as.data.frame(ode(Init, Time, feeding.markov.model, parms=params))

names(out) = c("time", "UF", "UP", "UN", "deltaF", "deltaP", "deltaN",
               "TF", "TP", "TN", "Dead")

out = out %>%
  mutate(prop_alive = 1-Dead)

plot(out[,2]/out[,12]~ out[,1], type="l", col="blue",lty="dashed",
     ylim=c(0,1), main="Compartment-specific death, 40% coverage",
     xlab="Time", ylab="Proportion of living mosquitoes")
lines(out[,3]/out[,12]~ out[,1], type="l", col="red", lty="dashed")
lines(out[,4]/out[,12]~ out[,1], type="l", col="green", lty="dashed")
lines(out[,8]/out[,12]~ out[,1], type="l", col="blue")
lines(out[,9]/out[,12]~ out[,1], type="l", col="red")
lines(out[,10]/out[,12]~ out[,1], type="l", col="green")

legend("topright", c("Untreated House", "Treated House", "Unfed",
                     "Part-fed",
                  "Fully fed"),
       lty=c(2, 1, 1, 1, 1),
       col=c("black", "black", "green", "red", "blue", "purple"))

