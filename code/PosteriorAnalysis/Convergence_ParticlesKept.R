library(Bolstad2)

setwd("~/Dropbox/DengueEntomologicalEffects/results/Posteriors")

# check at 25% of particles

load("SMC_Taper.RData")
samples1 = samples


load("SMC_Taper2.RData")
samples2 = samples

nrows = min(nrow(samples1), nrow(samples2))
samples1=head(samples1, nrows)
samples2=head(samples2, nrows)

conv = rep(NA, ncol(samples1))
for(i in 1:ncol(samples1)){
  
  theta = cbind(as.matrix(samples1[,i]), as.matrix(samples2[,i]))
  res = GelmanRubin(theta)
  conv[i] = res$R
}

mean(conv > 1.2)

# converged at 50 50, 50 25, 
# 75 25, 75 50

# at level 1.1
# .0227 50/25

par(mfrow=c(2,1))
hist(samples1$catch_prop, xlim=c(0, .06))
hist(samples2$catch_prop, xlim=c(0, .06))


##### combine the two samples into a single chain weighted by likelihood ####


samples = rbind(head(samples1, 100), head(samples2, 100))
lik.vec = rep(NA, 20)
for(i in 1:nrow(samples)){
  print(i)
  lik.vec[i] = likelihood(as.numeric(samples[i,]))
}

exp(lik.vec - max(lik.vec))
