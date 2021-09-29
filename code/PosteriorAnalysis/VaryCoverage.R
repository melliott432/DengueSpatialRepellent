#### load samples, prior, functions ####

load("~/Dropbox/DengueEntomologicalEffects/results/Posteriors/SMC_Taper4.RData")

setwd("~/Dropbox/DengueEntomologicalEffects/code")
source("functions.R")

#### make a data frame to hold all the information ####
coverage.vec = seq(0, 1, by=.05)
empty = rep(NA, length(coverage.vec))
coverage.df = data.frame(coverage=empty, 
                         modified_biting=empty,
                         modified_biting_lower=empty,
                         modified_biting_upper=empty,
                         modified_death=empty, 
                         modified_death_lower=empty, 
                         modified_death_upper=empty, 
                         time_untreated = empty,
                         time_untreated_lower = empty,
                         time_untreated_upper = empty,
                         time_treated=empty,
                         time_treated_lower=empty,
                         time_treated_upper=empty)

#### fill the data frame ####

for(i in 1:length(coverage.vec)){
  cov = coverage.vec[i]
  print(cov)
  
  coverage.df[i, "coverage"] = cov
  
  ac_calc_vec = rep(NA, nrow(samples))
  for(j in 1:nrow(samples)){
    ac_calc_vec[j] = with(samples[j,], ac_calc(af, ap, rho, C=cov, alpha, qt, qu, tau))
  }
  if(cov== 0){
    ac_calc_control=ac_calc_vec
  }
  relative_biting = ac_calc_vec / ac_calc_control
  coverage.df[i, "modified_biting"] = median(relative_biting)
  coverage.df[i, "modified_biting_lower"] = quantile(relative_biting, .025)
  coverage.df[i, "modified_biting_upper"] = quantile(relative_biting, .975)

  #### mosquito locations to get modified death rate
  mosquito.locations.treatment = matrix(rep(NA, nrow(samples)*9), nrow=nrow(samples),
                                        ncol=9)
  for(j in 1:nrow(samples)){
    mosquito.locations.treatment[j,] = with(samples[j,],  location_probability(af, ap,
                                                                               qu, tau,
                                                                               d, rho,
                                                                               qt, alpha,
                                                                               C=cov))
  } 
  samples$pi.U.treated =  mosquito.locations.treatment[,1] +
    mosquito.locations.treatment[,2]+mosquito.locations.treatment[,3]
  samples$pi.tau.treated =  mosquito.locations.treatment[,4] +
    mosquito.locations.treatment[,5]+mosquito.locations.treatment[,6]
  samples$pi.T.treated =  mosquito.locations.treatment[,7] +
    mosquito.locations.treatment[,8]+mosquito.locations.treatment[,9]

  #### modified death rate  
  
  gc_calc_vec = rep(NA, nrow(samples))
  for(j in 1:nrow(samples)){
    gc_calc_vec[j] = with(samples[j,], gc_calc(mu, qu, qt, gu, gt, gtau, tau, rho, C=cov,
                                               pi.tau.treated, pi.U.treated, pi.T.treated))
  }
  if(cov == 0){
    gc_calc_control=gc_calc_vec
  }
  relative_death = gc_calc_vec / gc_calc_control
  coverage.df[i, "modified_death"] = median(relative_death)
  coverage.df[i, "modified_death_lower"] = quantile(relative_death, .025)
  coverage.df[i, "modified_death_upper"] = quantile(relative_death, .975)
  
  
}

#### visualize ####
par(mfrow=c(1,2))
plot(modified_biting~coverage, coverage.df, ylim=c(.6, 1))
for(i in 1:nrow(coverage.df)){
  segments(coverage.df$coverage[i], coverage.df$modified_biting_lower[i],
           coverage.df$coverage[i],
            coverage.df$modified_biting_upper[i])
}
plot(modified_death~coverage, coverage.df, ylim=c(1,2.1))
for(i in 1:nrow(coverage.df)){
  segments(coverage.df$coverage[i], coverage.df$modified_death_lower[i],
           coverage.df$coverage[i],
           coverage.df$modified_death_upper[i])
}
