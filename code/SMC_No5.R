#### Performs Sequential Monte Carlo ####

rm(list = ls())
library(dplyr)
library(doParallel)


#### load files and functions ####

library(mvtnorm)
library(truncnorm)
source("functions.R")


#### load data from clinical trial ####

# import data: baseline
incidence.baseline <- read.csv("baseline.epi.csv")
hlc.baseline <- read.csv("baseline.hlc.csv")
parity.baseline = read.csv("baseline.parity.csv")
spor.rate.baseline = read.csv("baseline.spor.csv")
human.location <- read.csv("HumanLocation.csv")
cluster.info <- read.csv("cluster.info.csv")
coverage = read.csv("coverage.csv")

# standardize human location clusters
human.location = human.location %>%
  left_join(cluster.info, by = c("Cluster" = "study.cluster")) %>%
  rename(study.cluster = Cluster) %>%
  rename(Cluster = analysis.cluster) %>%
  dplyr::select(-c("X", "study.cluster"))


# add human location data
hlc.baseline = hlc.baseline %>%
  left_join(human.location, by = "Cluster")


# intervention
incidence.intervention <- read.csv("intervention.epi.csv")
hlc.intervention <- read.csv("intervention.hlc.csv")
parity.intervention = read.csv("intervention.parity.csv")
spor.rate.intervention = read.csv("intervention.spor.csv")


hlc.intervention = hlc.intervention %>%
  left_join(human.location, by = "Cluster")


# indicate t and c clusters
treatment_clusters = c(1, 5, 6, 7, 8, 9)
control_clusters = c(2,3,4, 10, 11, 12)

# import prior samples
load("PriorGeneratorOutputInterim.RData")

prior.samps = prior.samples

# remove any NAs
prior.samps = prior.samps[rowSums(is.na(prior.samps)) == 0,]


# add truncated normal priors for multiplicative values
prior.samps$a.mult = rtruncnorm(nrow(prior.samps), a = 0,
                                b = Inf,
                                mean = 1,
                                sd = .2)

prior.samps$g.mult = rtruncnorm(nrow(prior.samps), a = 0,
                                b = Inf,
                                mean = 1,
                                sd = .2)




# add a instead of f and q
prior.samps = prior.samps %>%
  mutate(a=f*q) %>%
  dplyr::select(-c( f, q))

# arrange
prior.samps = prior.samps %>%
  dplyr::select(c(paste0("lambda", 1:12), au=a, a.mult, gu=g, g.mult, b,
                  paste0("c", 1:12), n, r))

# set q 
q=.25


#### define the log likelihood ####
likelihood = function(params){
  
  # separate out the parameters and name them for easier referencing
  lambda = params[1:12]
  
  au = params[13]
  a.mult = params[14]
  at = au*a.mult
  
  gu = params[15]
  g.mult = params[16]
  gt = gu*g.mult
  
  b = params[17]
  c = params[18:29]
  n = params[30]
  r = params[31]
  
  # coverage for each cluster
  psi = coverage$coverage[treatment_clusters]
  a_psi = au*(1-psi) + at*psi
  
  # calcualte the parity rate for all at baseline and for control clusters during intervention
  parity.rate.baseline = parity.rate.intervention = rep(parity_control(au,gu,q), 12)
  # update intervention parity for treated clusters
  parity.rate.intervention[treatment_clusters] = parity_treatment(au, a.mult,
                                                                  gu, g.mult, q, psi)
  
  # same for incidence, but with different lambda and c values- gives 12 different values of incidence
  inc.baseline = inc.intervention = incidence_control(au, gu, lambda, b, c, n, r)
  incidence_treated_clusters = rep(NA, length(treatment_clusters))
  for(i in 1:length(treatment_clusters)){
    incidence_treated_clusters[i] = incidence_treated(au, psi[i],
                                                      lambda[treatment_clusters[i]],
                                                      b, c[treatment_clusters[i]], gu,
                                                      g.mult, n, r, a.mult)
  }
  
  inc.intervention[treatment_clusters] = incidence_treated_clusters
  
  # daily Human biting rate as calculated from proposed values of parameters
  hbr.baseline = hbr.intervention = human_biting_rate_control(au, gu, lambda)
  hbr.intervention[treatment_clusters] = human_biting_rate_treated_intervention(au, a.mult, psi,
                                                                                gu, g.mult, 
                                                                                lambda[treatment_clusters])
  
  # sporozoite rate, similar to incidence
  sporo.rate.baseline = sporo.rate.intervention = prevalence_mosquitoes_control(au, gu, lambda,
                                                                                b, c, n, r)
  
  
  mosquito_prev_treated_clusters = rep(NA, length(treatment_clusters))
  for(i in 1:length(treatment_clusters)){
    mosquito_prev_treated_clusters[i] =  prevalence_mosquitoes_treatment(au, psi[i],a.mult, c[treatment_clusters[i]], 
                                                                         gu, g.mult, n, 
                                                                         lambda[treatment_clusters[i]], r, b)
    
    
  }
  
  sporo.rate.intervention[treatment_clusters] = mosquito_prev_treated_clusters
  
  # daily indoor biting rate = proportion of people inside * number of bites per hour
  indoor.biting.rate.baseline = human.location$prop.inside * hbr.baseline
  indoor.biting.rate.intervention = human.location$prop.inside * hbr.intervention
  
  # daily outdoor biting rate, calculation the same
  outdoor.biting.rate.baseline = human.location$prop.outside * hbr.baseline
  outdoor.biting.rate.intervention = human.location$prop.outside * hbr.intervention
  
  # likelihood for parity
  parity.likelihood.baseline = sum(dbinom(x = parity.baseline$number.parous,
                                          size = parity.baseline$number.tested,
                                          prob = parity.rate.baseline, log = T)[-5])
  parity.likelihood.intervention = sum(dbinom(x = parity.intervention$number.parous,
                                              size = parity.intervention$number.tested,
                                              prob = parity.rate.intervention, log = T)[-5])
  
  # likelihood for hlc
  hlc.baseline = hlc.baseline %>%
    filter(Cluster != 5)
  hlc.intervention = hlc.intervention %>%
    filter(Cluster != 5)
  indoor.hlc.likelihood.baseline = sum(dpois(hlc.baseline$inside,
                                             lambda = indoor.biting.rate.baseline[hlc.baseline$Cluster]*
                                               # multiply by 5/6 because HLC was only performed for 50 minutes each hour\
                                               # so the biting rate is only 5/6 of what it would be.
                                               # did this in lambda
                                               #  to avoid non-integer in the Poisson likleihood
                                               5/6, log = T))
  
  indoor.hlc.likelihood.intervention = sum(dpois(hlc.intervention$inside,
                                                 indoor.biting.rate.intervention[hlc.intervention$Cluster]*
                                                   5/6, log = T))
  
  outdoor.hlc.likelihood.baseline = sum(dpois(hlc.baseline$outside,
                                              outdoor.biting.rate.baseline[hlc.baseline$Cluster]*5/6,
                                              log = T))
  
  outdoor.hlc.likelihood.intervention = sum(dpois(hlc.intervention$outside,
                                                  outdoor.biting.rate.intervention[hlc.intervention$Cluster]*
                                                    5/6, log = T))
  
  # likelihood for incidence
  incidence.baseline = incidence.baseline %>%
    filter(Cluster != 5)
  incidence.intervention = incidence.intervention %>%
    filter(Cluster != 5)
  incidence.likelihood.baseline = sum(dpois(x = incidence.baseline$total.infections,
                                            # total infections per person over entire arm
                                            lambda = (inc.baseline[incidence.baseline$Cluster] +
                                                        # add a very small number because dpois would return 0 if 
                                                        # not logged and -Inf when logged
                                                        .0000000000000000001) *
                                              # divide folloup time by 365 because it is in days and incidence is yearly
                                              incidence.baseline$followup.time / 365, log = T))
  incidence.likelihood.intervention = sum(dpois(incidence.intervention$total.infections,
                                                (inc.intervention[incidence.intervention$Cluster] +
                                                   .0000000000000000001) *
                                                  incidence.intervention$followup.time / 365,
                                                log = T))
  
  # spor rate likelihood
  spor.rate.likelihood.baseline = sum(dbinom(spor.rate.baseline$total.positive,
                                             size = spor.rate.baseline$total.caught,
                                             # add a small number so prob is not 0
                                             prob = sporo.rate.baseline + .00000000000000000000000001,
                                             log = T)[-5])
  
  spor.rate.likelihood.intervention= sum(dbinom(spor.rate.intervention$total.positive,
                                                size = spor.rate.intervention$total.caught,
                                                prob = sporo.rate.intervention + 
                                                  .00000000000000000000000001, log = T)[-5])
  
  return(
    # comment out sporozoite rate if needed
    parity.likelihood.baseline + 
      indoor.hlc.likelihood.baseline + 
      outdoor.hlc.likelihood.baseline + 
      incidence.likelihood.baseline +
      spor.rate.likelihood.baseline + 
      parity.likelihood.intervention + 
      indoor.hlc.likelihood.intervention +  
      outdoor.hlc.likelihood.intervention + 
      incidence.likelihood.intervention  +  
      spor.rate.likelihood.intervention
  )
  
}


#### Function to perform sequential Monte-Carlo ####
samples = prior.samps
sequential_mc_one_round = function(samples, prop.particles=.1, target.samples=.9){
  
  target.particles = round(prop.particles * nrow(samples))
  target.samples = target.samples
  
  
  # log.ll = rep(NA, nrow(samples))
  # for(i in 1:nrow(samples)){
  #   print(i)
  #   log.ll[i] = likelihood(as.numeric(samples[i,]))
  # }
  
  log.ll = foreach(i=1:nrow(samples), .combine='c', .export=ls(globalenv()),
                   .packages = c("dplyr","mvtnorm","doParallel")) %dopar% {
                     likelihood(as.numeric(samples[i,]))
                   }
  
  
  samples$likelihood = log.ll
  
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
  
  
  # If there is an NA in prob, set it to 0
  # wts[is.na(wts)] = 0
  
  ind = sample(1:nrow(samples),
         size = nrow(samples), replace = T,
         prob = wts)
  resampled = samples[ind,-32]
  
  prop_kept = length(unique(ind)) / nrow(samples)
  
  
  
  mu = colMeans(resampled)
  sigma = var(resampled)
  
  
  new_df = data.frame(rmvnorm(100000, mu, sigma))
  
  # filter out any rows with any negative values or 
  # values outside (0,1) for probabilities
  new_df = new_df[rowSums(new_df < 0) == 0,]
  less_than_1_cols = c(15,17, 18:28)
  new_df = new_df[rowSums(new_df[, less_than_1_cols] > 1) == 0, ]
  
  return(list(resampled, mu, sigma, const, new_df, prop_kept))
}

#### Perform 1000 times ####

cl = makeCluster(20)
registerDoParallel(cl)

times = 100*99
mu.list = vector(mode = "list", length=times)
sigma.list = vector(mode = "list", length=times)
prop_kept_vec = vector(mode="list", length=times)

const.vec = rep(NA, times)

# Inititalize samples 
samples = prior.samps
current.times = 1
while(current.times <= times){
  prop.particles = 1-ceiling(current.times/100)/100
  res = sequential_mc_one_round(samples, prop.particles)
  samples = res[[5]]
  mu.list[[current.times]] = res[[2]]
  sigma.list[[current.times]] = res[[3]]
  const.vec[current.times] = res[[4]]
  prop_kept_vec[[current.times]] = res[[6]]
  current.times = current.times + 1
  save(samples, mu.list, sigma.list, const.vec,
       file = "SMC_No5.RData")
  samples = res[[1]]
}

stopCluster(cl)