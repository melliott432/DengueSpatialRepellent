#### Performs Sequential Monte Carlo ####
rm(list=ls())

library(dplyr)
library(mvtnorm)
library(doParallel)


#### load files and functions ####

source("functions.R")
source("LikelihoodFunction.R")

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

#### Priors ####

# import prior samples
load("PriorSamples.RData")


# arrange so each row is in the right order for the likelihood function
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
                  paste0("X_", 1:26)), gtau,
                phi, d, catch_prop)


# set coverage
C = 0.45
p = 5.2



#### Function to perform sequential Monte-Carlo ####

sequential_mc_one_round = function(samples){
  
  target.particles = round(0.75 * nrow(samples))
  target.samples = 0.25
  
  
  
  
  # log.ll = rep(NA, nrow(samples))
  # 
  # for(i in 1:nrow(samples)){
  # log.ll[i] = likelihood(as.numeric(samples[i,]))
  #  }
  
  log.ll = foreach(i=1:nrow(samples), .combine='c', .export=ls(globalenv()),
                   .packages = c("dplyr","mvtnorm","doParallel", "deSolve",
                                 "msm", "markovchain")) %dopar% {
                                   likelihood(as.numeric(samples[i,]))
                                 }
  samples$likelihood = log.ll
  # a higher constant penalizes lower likelihoods more
  find.const = function(const){
    wts = exp(
      const *
        (samples$likelihood - max(samples$likelihood)))
    wts = wts / sum(wts)
    (cumsum(sort(wts,decreasing=T))[target.particles] - target.samples) ^ 2
  }
  
  const = optimize(
    find.const,
    interval=c(1e-10,1e-3),
    tol=.Machine$double.eps^0.5)$minimum
  
  wts = exp(
    const * (
      samples$likelihood -
        max(samples$likelihood)))
  wts = wts / sum(wts)
  
  wts[is.na(wts)] = 0
  resampled = samples[sample(1:nrow(samples),
                             size = nrow(samples), replace = T,
                             prob = wts),-45]
  
  
  
  mu = colMeans(resampled)
  sigma = var(resampled)
  
  
  new_df = data.frame(rmvnorm(100000, mu, sigma))
  
  # filter out any rows with any negative values or 
  # values outside (0,1) for probabilities
  new_df = new_df[rowSums(new_df < 0) == 0,]
  less_than_1_cols = c(3, 4, 5, 10, 11, 14, 15:40, 41,44)
  new_df = new_df[rowSums(new_df[, less_than_1_cols] > 1) == 0, ]
  
  return(list(new_df, mu, sigma, const))
}

#### register parallel ####

cl = makeCluster(20)
registerDoParallel(cl)



#### Perform 100 times ####
times = 100000
mu.list = vector(mode = "list", length=times)
sigma.list = vector(mode = "list", length=times)

const.vec = rep(NA, times)

# Inititalize samples 
samples = prior.samps
current.times = 1
while(current.times <= times){
  res = sequential_mc_one_round(samples)
  samples = res[[1]]
  mu.list[[current.times]] = res[[2]]
  sigma.list[[current.times]] = res[[3]]
  const.vec[current.times] = res[[4]]
  current.times = current.times + 1
  save(samples, mu.list, sigma.list, const.vec,
       file = "SMC_75_25.RData")
}

stopCluster(cl)

save(samples, mu.list, sigma.list, const.vec,
     file = "SMC_75_25.RData")