library(Bolstad2)

#### import the posteriors ####
rm(list=ls())
setwd("~/Dropbox/DengueEntomologicalEffects/results/Posteriors")

load("SequentialMonteCarlo.RData")
samples1=samples
load("SequentialMonteCarlo2.RData")
samples2=samples

number.samples = min(nrow(samples1), nrow(samples2))
theta = cbind(as.matrix(head(samples1$r, number.samples)), as.matrix(head(samples2$r, number.samples)))
res = GelmanRubin(theta)
res



samples1=samples
mu.list1=mu.list

for(i in 2:10){
  print(i)
  file_name = paste0("SequentialMonteCarlo", i, ".RData")
  load(file_name)
  load("SequentialMonteCarlo2.RData")
  assign(paste0("samples", i), samples)
  assign(paste0("mu.list", i), mu.list)
}

#### Mu list ####

# data frame to hold mu over time: try lambda, a.mult, g.mult, and r

n.df = matrix(rep(NA, 10*1700), nrow=1700, ncol=10)
rho.df = matrix(rep(NA, 10*1700), nrow=1700, ncol=10)
X_12.df = matrix(rep(NA, 10*1700), nrow=1700, ncol=10)
qu.df = matrix(rep(NA, 10*1700), nrow=1700, ncol=10)

for(i in 1:1700){
  for(j in 1:10){
    n.df[i,j] = get(paste0("mu.list", j))[[i]]["n"]
    rho.df[i,j] = get(paste0("mu.list", j))[[i]]["rho"]
    X_12.df[i,j] = get(paste0("mu.list", j))[[i]]["X_12"]
   qu.df[i,j] = get(paste0("mu.list", j))[[i]]["qu"]
  }
}


png("~/Dropbox/DengueEntomologicalEffects/results/ManuscriptPlots/Convergence.png")
par(mfrow=c(2,2))

plot(n.df[,1], col="red", type="l", main="Convergence of n")
lines(n.df[,2], col="blue", type="l")
lines(n.df[,3], col="yellow", type="l")
lines(n.df[,4], col="green", type="l")
lines(n.df[,5], col="orange", type="l")
lines(n.df[,6], col="black", type="l")
lines(n.df[,7], col="brown", type="l")
lines(n.df[,8], col="purple", type="l")
lines(n.df[,9], col="pink", type="l")
lines(n.df[,10], col="darkred", type="l")


plot(rho.df[,1], col="red", type="l", main="Convergence of repellency probability")
lines(rho.df[,2], col="blue", type="l")
lines(rho.df[,3], col="yellow", type="l")
lines(rho.df[,4], col="green", type="l")
lines(rho.df[,5], col="orange", type="l")
lines(rho.df[,6], col="black", type="l")
lines(rho.df[,7], col="brown", type="l")
lines(rho.df[,8], col="purple", type="l")
lines(rho.df[,9], col="pink", type="l")
lines(rho.df[,10], col="darkred", type="l")

plot(X_12.df[,1], col="red", type="l", main="Convergence of prevalence Cluster 12", ylim=c(0,.2))
lines(X_12.df[,2], col="blue", type="l")
lines(X_12.df[,3], col="yellow", type="l")
lines(X_12.df[,4], col="green", type="l")
lines(X_12.df[,5], col="orange", type="l")
lines(X_12.df[,6], col="black", type="l")
lines(X_12.df[,7], col="brown", type="l")
lines(X_12.df[,8], col="purple", type="l")
lines(X_12.df[,9], col="pink", type="l")
lines(X_12.df[,10], col="darkred", type="l")

plot(qu.df[,1], col="red", type="l", main="Convergence of qu")
lines(qu.df[,2], col="blue", type="l")
lines(qu.df[,3], col="yellow", type="l")
lines(qu.df[,4], col="green", type="l")
lines(qu.df[,5], col="orange", type="l")
lines(qu.df[,6], col="black", type="l")
lines(qu.df[,7], col="brown", type="l")
lines(qu.df[,8], col="purple", type="l")
lines(qu.df[,9], col="pink", type="l")
lines(qu.df[,10], col="darkred", type="l")

dev.off()

library(Bolstad2)


theta = cbind(as.matrix(samples2$r), as.matrix(samples3$r))
res = GelmanRubin(theta)
res
