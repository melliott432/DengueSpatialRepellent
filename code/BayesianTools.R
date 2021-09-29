rm(list = ls())
library(dplyr)
library(doParallel)
library(BayesianTools)


library(mvtnorm)
library(truncnorm)
source("functions.R")


#### load data from clinical trial ####

epi = read.csv(file="epi.csv")
parity_baseline = read.csv(file="parity_baseline.csv")
parity_intervention = read.csv(file="parity_intervention.csv")
load("treatedclusters.RData")
abundance_baseline = read.csv(file="abundance_baseline.csv")
abundance_intervention_control = read.csv(file="abundance_intervention_control.csv")
abundance_intervention_treatment = read.csv(file="abundance_intervention_treatment.csv")
bloodmeal_baseline = read.csv(file="bloodmeal_baseline.csv")
bloodmeal_intervention = read.csv(file="bloodmeal_intervention.csv")

load("PriorSamples.RData")


prior.samps = prior.samps %>%
  dplyr::select(c(alpha,
                  qt,
                  gt,
                  mu,
                  rho,
                  lambda,
                  af,
                  ap,
                  qu ,
                  c,
                  gu,
                  n,
                  tau,
                  
                  # human parameters
                  b,
                  r,
                  paste0("X_", 1:26)), gtau,
                phi, d, catch_prop)


# set coverage
C = 0.45
p = 5.2

source("LikelihoodFunction.R")

#### fit a multivariate normal to prior distribution ####

mu = colMeans(prior.samps)
sigma = var(prior.samps)

priorDensity = function(params){
  return(dmvnorm(params, mu, sigma, log=T))
}

sampler = function(n=1){
  value_chosen=F
  while(value_chosen==F){
    potential.value = rmvnorm(n, mu, sigma)
    if(sum(potential.value < 0) == 0){
      if(sum(potential.value[c(3, 4, 5, 10, 11, 14, 16:41, 42, 45)]>1) == 0){
        value_chosen=T
        return(potential.value)
      }
    }
  }
}


bayesianPrior = createPrior(priorDensity,sampler, lower=rep(0, 45), upper=c(
  2, 1000, 1, 1, 1, 1000, 1000, 1000, 1000, 1, 1, 100, 1000, 1, 1000, rep(1, 26),
  1, 1000, 1000, 1
))


BayesianSetup = createBayesianSetup(likelihood=likelihood, prior=bayesianPrior)
out=runMCMC(BayesianSetup, settings=list(iterations=100000))
codaObject = data.frame(getSample(out, start = 500, thin=10))
names(codaObject) = names(prior.samps)

save(codaObject, out, file="BayesianTools.RData")
