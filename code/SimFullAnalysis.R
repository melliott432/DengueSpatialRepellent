library(dplyr)

# load the data sets and parameter values used to make the simulation
setwd("~/OneDrive/Documents/ND/Research/DengueSpatialRepellent/simulations_data")
load("SimFull0601.RData")
load("20SimulatedDatasets.RData")

#### Comparison of posterior to true value ####
npars = 8

# reorder combos to try
combos_to_try = head(combos_to_try, 20)
names(combos_to_try) = sub("_true", "", names(combos_to_try))

combos_to_try = combos_to_try %>%
  dplyr::select(lambda, au, a.mult, g.mult, q.mult, n, X, gu)

# determine what proportion of the posterior contained the true value for each parameter
prop_true_value = rep(0,npars)
for(par in 1:npars){
  for(i in 1:20){
    if (quantile(result[[i]][[1]][,par], .005) < combos_to_try[i, par] & 
        quantile(result[[i]][[1]][,par], .995) > combos_to_try[i,par]){
      prop_true_value[par] = prop_true_value[par]+1
    }
  }
}

prop_true_value = prop_true_value / 20

setwd("~/OneDrive/Documents/ND/Research/SR/plots")
png("2022-06-01-Simulation_ValueComparison_7.png", height=1000, width=1000)
# plot posterior medians vs true parameter
par(mfrow=c(2,4), cex=1.1, mar=c(3, 3, 3, 1))
lows = rep(0, npars)
tops = c(100, 
         5, 
         2, 
         2, 
         2, 
         #10,10,10, 
         20, 
         #2, 
         1, 
         1)
         #1, 1, 1, 1)
param_names = c("Mosquito Emergence", 
                "Baseline biting rate", 
                "Multiplier on biting",
                "Multiplier on mortality",
                "Multiplier on Exit Rate",
                #"Dispersion",
                #"Baseline Exit Rate", 
                #"Entrance Rate", 
                "EIP",
                #"Digestion Rate",
                "Arbovirus Prevalence",
                "Baseline Mortality")
                #"Mosquito to Human Transmission",
                #"Human to Mosquito Transmission", 
                #"Repellency",
                #"Proportion Mosquitoes Caught")
for(par in 1:npars){
  plot(0, type="n", xlim=c(lows[par], tops[par]), ylim=c(lows[par], tops[par]), xlab="True",
       ylab="Inferred",
       main=param_names[par])
  for(i in 1:20){
    points(combos_to_try[i,par], median(result[[i]][[1]][,par]), col="red", pch=19)
  }
  abline(0, 1)
  # add bars for posterior distribution
  for(i in 1:20){
    segments(x0=combos_to_try[i,par],x1=combos_to_try[i,par], 
             y0 = quantile(result[[i]][[1]][,par], .025),
             y1 = quantile(result[[i]][[1]][,par], .975))
  }
}
dev.off()


setwd("~/OneDrive/Documents/ND/Research/DengueSpatialRepellent/code")
source("functions.R")

#### load trial data ####

setwd("~/OneDrive/Documents/ND/Research/DengueSpatialRepellent/likelihood_data")
bloodmeal_intervention = read.csv("bloodmeal_intervention.csv")
bloodmeal_baseline = read.csv("bloodmeal_baseline.csv")

# take total full to be the sum of hald and full
bloodmeal_intervention = bloodmeal_intervention %>%
  mutate(total_fed = total_full + total_half) %>%
  dplyr::select(-c(total_full, total_half))

bloodmeal_baseline = bloodmeal_baseline %>%
  mutate(total_fed = total_full + total_half) %>%
  dplyr::select(-c(total_full, total_half))

abundance_baseline = read.csv("abundance_baseline.csv")
abundance_intervention_control = read.csv("abundance_intervention_control.csv")
abundance_intervention_treatment = read.csv("abundance_intervention_treatment.csv")

parity_baseline = read.csv("parity_baseline.csv")
parity_intervention = read.csv("parity_intervention.csv")

epi = read.csv("epi.csv")

load("treatedclusters.RData")
C=.5
p=5.2

#### Calculate necessary quantities ####

# proportional locations calculated from combos
posterior_list = list()

for(i in 1:20){
  posterior = result[[i]][[1]]
  posterior = head(posterior, 100)
  print(i)
  qu = posterior$qu
  gu = posterior$gu
  tau = posterior$tau
  d = posterior$d
  b = posterior$b
  c = posterior$c
  n = posterior$n
  X = posterior$X
  rho = posterior$rho
  q.mult = posterior$q.mult
  a.mult = posterior$a.mult
  g.mult = posterior$g.mult
  locations_control = list()
  locations_treatment = list()
  for(j in 1:nrow(posterior)){
    locations_control[[j]] = with(posterior[j,], location_probability(posterior$au[j], qu,
                                                                      tau, d,
                                                                      rho, q.mult,
                                                                      a.mult, C=0))
    locations_treatment[[j]] = with(posterior[j,], location_probability(posterior$au[j],
                                                                        qu, tau, d,
                                                                        rho, q.mult,
                                                                        a.mult,
                                                                        C))
  }
  # get posterior distributions of mosquito locations
  posterior$UF_control = unlist(lapply(locations_control, function(l) l[1]))
  posterior$UN_control = unlist(lapply(locations_control, function(l) l[2]))
  posterior$DF_control = unlist(lapply(locations_control, function(l) l[3]))
  posterior$DN_control = unlist(lapply(locations_control, function(l) l[4]))
  
  posterior$UF_treatment = unlist(lapply(locations_treatment, function(l) l[1]))
  posterior$UN_treatment = unlist(lapply(locations_treatment, function(l) l[2]))
  posterior$DF_treatment = unlist(lapply(locations_treatment, function(l) l[3]))
  posterior$DN_treatment = unlist(lapply(locations_treatment, function(l) l[4]))
  posterior$TF_treatment = unlist(lapply(locations_treatment, function(l) l[5]))
  posterior$TN_treatment = unlist(lapply(locations_treatment, function(l) l[6]))
  
  # get posterior distributions of control and treatment modified biting rates
  posterior$overall_death_control = gc_calc(gu, g.mult,
                                            pi.tau=posterior$DN_control+posterior$DF_control,
                                            pi.U=posterior$UN_control+posterior$UF_control,
                                            pi.T=0)
  posterior$overall_death_treatment = gc_calc(gu, g.mult,
                                              posterior$DN_treatment+posterior$DF_treatment,
                                              posterior$UN_treatment+posterior$UF_treatment,
                                              pi.T=posterior$TN_treatment+
                                                posterior$TF_treatment)
  
  posterior$overall_biting_control = ac_calc(posterior$au, a.mult, d,
                                             posterior$UN_control,
                                             0,
                                             posterior$DN_control)
  posterior$overall_biting_treatment = ac_calc(posterior$au, a.mult, d,
                                               posterior$UN_treatment,
                                               posterior$TN_treatment,
                                               posterior$DN_treatment)
  
  # find a distribution over parity
  posterior$parity_control = parity(posterior$overall_biting_control,
                                    posterior$overall_death_control)
  posterior$parity_treatment = parity(posterior$overall_biting_treatment,
                                      posterior$overall_death_treatment)
  
  # find a distribution over force of infection
  posterior$foi_control = force_of_infection(b, 
                                             posterior$lambda,
                                             posterior$overall_biting_control,
                                             c, X, 
                                             posterior$overall_death_control, n)
  posterior$foi_treatment = force_of_infection(b, posterior$lambda,
                                               posterior$overall_biting_treatment,
                                               c, X, posterior$overall_death_treatment, n)
  
  # find a distribution over expected number of mosquitoes in a house
  posterior$m_control = posterior$lambda / posterior$overall_death_control
  posterior$m_treatment = posterior$lambda / posterior$overall_death_treatment
  posterior$f_control = posterior$m_control * p * (posterior$UF_control+posterior$UN_control)
  posterior$f_treatment = posterior$m_treatment * p * (posterior$TF_treatment+
                                                         posterior$TN_treatment)
  
  # find a distribution over the proportion of mosquitoes that are blood fed
  posterior$proportion_fed_control = posterior$UF_control / (posterior$UF_control + 
                                                               posterior$UN_control)
  posterior$proportion_fed_treatment = posterior$TF_treatment / (posterior$TF_treatment + 
                                                                   posterior$TN_treatment)
  
  result[[i]][[1]] = posterior
}

#### PPD for parity ####
par_control_post_med = rep(NA, length(result))
par_control_post_lower = rep(NA, length(result))
par_control_post_upper = rep(NA, length(result))

par_int_control_post_med = rep(NA, length(result))
par_int_control_post_lower = rep(NA, length(result))
par_int_control_post_upper = rep(NA, length(result))

par_int_treatment_post_med = rep(NA, length(result))
par_int_treatment_post_lower = rep(NA, length(result))
par_int_treatment_post_upper = rep(NA, length(result))

for(i in 1:length(result)){
  posterior = result[[i]][[1]]
  
  # use each of these to predict a new data point
  post_preds_baseline = rbinom(1000, size=sum(parity_baseline$total_caught), 
                               prob=posterior$parity_control)
  post_preds_intervention_control = rbinom(1000,
                                           size=sum(parity_intervention[parity_intervention$cluster_trt == "C",]$total_caught), 
                                           prob=posterior$parity_control)
  post_preds_intervention_treatment = rbinom(1000,
                                             size=sum(parity_intervention[
                                               parity_intervention$cluster_trt == "T",
                                             ]$total_caught), 
                                             prob=posterior$parity_treatment)
  
  # add to vector
  par_control_post_med[i] = median(post_preds_baseline) 
  par_control_post_lower[i] =  quantile(post_preds_baseline, .025) 
  par_control_post_upper[i] =  quantile(post_preds_baseline, .975) 
  
  
  par_int_control_post_med[i] = median(post_preds_intervention_control) 
  par_int_control_post_lower[i] =  quantile(post_preds_intervention_control, .025) 
  par_int_control_post_upper[i] =  quantile(post_preds_intervention_control, .975) 
  
  par_int_treatment_post_med[i] = median(post_preds_intervention_treatment) 
  par_int_treatment_post_lower[i] =  quantile(post_preds_intervention_treatment, .025) 
  par_int_treatment_post_upper[i] =  quantile(post_preds_intervention_treatment, .975) 
  
}

# plot vs the true values
true_par_control = rep(NA, length(result))
for(i in 1:length(result)){
  true_par_control[i] = simulated_datasets[[i]][[2]][[1]]
}


true_par_int_control = rep(NA, length(result))
for(i in 1:length(result)){
  true_par_int_control[i] = simulated_datasets[[i]][[2]][[2]]
}


true_par_int_treatment = rep(NA, length(result))
for(i in 1:length(result)){
  true_par_int_treatment[i] = simulated_datasets[[i]][[2]][[3]]
}

# plot the true values vs the actual
setwd("~/Dropbox/DengueEntomologicalEffects/Analyses/BayesianTools/results/plots")
png("PPD_parity_7.png", height=800, width=800)
par(mfrow=c(2,2), mar=c(5,5,4,3), cex=1.1)

# baseline
plot(0, type="n", xlim=c(150, 500), ylim=c(150, 500), xlab="Simulated Data", 
     ylab="Posterior Prediction", main="Baseline")
points(true_par_control, par_control_post_med)
abline(a=0, b=1)
segments(x0=true_par_control, y0=par_control_post_lower, x1=true_par_control,
         y1=par_control_post_upper)

# intervention control
plot(0, type="n", xlim=c(600,1600), ylim=c(600,1600), xlab="Simulated Data", 
     ylab="Posterior Prediction", main="Intervention, Control")
points(true_par_int_control, par_int_control_post_med)
abline(a=0, b=1)
segments(x0=true_par_int_control, y0=par_int_control_post_lower, x1=true_par_int_control,
         y1=par_int_control_post_upper)

# intervention treatment
plot(0, type="n", xlim=c(500, 1400), ylim=c(500, 1400), xlab="Simulated Data", 
     ylab="Posterior Prediction", main="Intervention, Treatment")
points(true_par_int_treatment, par_int_treatment_post_med)
abline(a=0, b=1)
segments(x0=true_par_int_treatment, y0=par_int_treatment_post_lower, x1=true_par_int_treatment,
         y1=par_int_treatment_post_upper)
dev.off()

#### PPD for bloodfeeding ####

blood_control_post_med = rep(NA, length(result))
blood_control_post_lower = rep(NA, length(result))
blood_control_post_upper = rep(NA, length(result))

blood_int_control_post_med = rep(NA, length(result))
blood_int_control_post_lower = rep(NA, length(result))
blood_int_control_post_upper = rep(NA, length(result))

blood_int_treatment_post_med = rep(NA, length(result))
blood_int_treatment_post_lower = rep(NA, length(result))
blood_int_treatment_post_upper = rep(NA, length(result))

for(i in 1:length(result)){
  posterior = result[[i]][[1]]
  
  # use each of these to predict a new data point
  post_preds_baseline = rbinom(1000, size=sum(bloodmeal_baseline$total_evaluated), 
                               prob=posterior$proportion_fed_control)
  post_preds_intervention_control= rbinom(1000, size=sum(bloodmeal_intervention[bloodmeal_intervention$cluster_trt == "C",]$total_evaluated), 
                                          prob=posterior$proportion_fed_control)
  post_preds_intervention_treatment = rbinom(1000, size=sum(bloodmeal_intervention[bloodmeal_intervention$cluster_trt == "T",]$total_evaluated), 
                                             prob=posterior$proportion_fed_treatment)
  
  # add to vector
  blood_control_post_med[i] = median(post_preds_baseline) 
  blood_control_post_lower[i] =  quantile(post_preds_baseline, .025) 
  blood_control_post_upper[i] =  quantile(post_preds_baseline, .975) 
  
  
  blood_int_control_post_med[i] = median(post_preds_intervention_control) 
  blood_int_control_post_lower[i] =  quantile(post_preds_intervention_control, .025) 
  blood_int_control_post_upper[i] =  quantile(post_preds_intervention_control, .975) 
  
  blood_int_treatment_post_med[i] = median(post_preds_intervention_treatment) 
  blood_int_treatment_post_lower[i] =  quantile(post_preds_intervention_treatment, .025) 
  blood_int_treatment_post_upper[i] =  quantile(post_preds_intervention_treatment, .975) 
  
}

# plot vs the true values
true_blood_control = rep(NA, length(result))
for(i in 1:length(result)){
  true_blood_control[i] = simulated_datasets[[i]][[4]][[1]]
}



true_blood_int_control = rep(NA, length(result))
for(i in 1:length(result)){
  true_blood_int_control[i] = simulated_datasets[[i]][[4]][[2]]
}


true_blood_int_treatment = rep(NA, length(result))
for(i in 1:length(result)){
  true_blood_int_treatment[i] = simulated_datasets[[i]][[4]][[3]]
}

# plot the true values vs the actual
setwd("~/Dropbox/DengueEntomologicalEffects/Analyses/BayesianTools/results/plots")
png("PPD_blood_7.png", height=800, width=800)
par(mfrow=c(2,2), mar=c(5,5,4,3), cex=1.1)

# baseline
plot(0, type="n", xlim=c(0, 3500), ylim=c(0, 3500), xlab="Simulated Data", 
     ylab="Posterior Prediction", main="Baseline")
points(true_blood_control, blood_control_post_med)
abline(a=0, b=1)
segments(x0=true_blood_control, y0=blood_control_post_lower, x1=true_blood_control,
         y1=blood_control_post_upper)

# intervention control
plot(0, type="n", xlim=c(0,13000), ylim=c(0, 13000), xlab="Simulated Data", 
     ylab="Posterior Prediction", main="Intervention, Control")
points(true_blood_int_control, blood_int_control_post_med)
abline(a=0, b=1)
segments(x0=true_blood_int_control, y0=blood_int_control_post_lower, x1=true_blood_int_control,
         y1=blood_int_control_post_upper)

# intervention treatment
plot(0, type="n", xlim=c(0, 13000), ylim=c(0 ,13000), xlab="Simulated Data", 
     ylab="Posterior Prediction", main="Intervention, Treatment")
points(true_blood_int_treatment, blood_int_treatment_post_med)
abline(a=0, b=1)
segments(x0=true_blood_int_treatment, y0=blood_int_treatment_post_lower, x1=true_blood_int_treatment,
         y1=blood_int_treatment_post_upper)
dev.off()

#### PPD for abundance ####

ab_control_post_med = rep(NA, length(result))
ab_control_post_lower = rep(NA, length(result))
ab_control_post_upper = rep(NA, length(result))

ab_int_control_post_med = rep(NA, length(result))
ab_int_control_post_lower = rep(NA, length(result))
ab_int_control_post_upper = rep(NA, length(result))

ab_int_treatment_post_med = rep(NA, length(result))
ab_int_treatment_post_lower = rep(NA, length(result))
ab_int_treatment_post_upper = rep(NA, length(result))

for(i in 1:length(result)){
  print(i)
  posterior = result[[i]][[1]]
  catch_prop = posterior$catch_prop
  phi = posterior$phi
  
  # use each of these to predict a new data point
  post_preds_baseline = unlist(lapply(1:nrow(posterior), function(row) {
    mean(rnbinom(nrow(abundance_baseline),mu = posterior$f_control[row]*
                   catch_prop, size=phi))}
  ))
  post_preds_intervention_control = unlist(lapply(1:nrow(posterior), function(row) {
    mean(rnbinom(nrow(abundance_intervention_control),
                 mu = posterior$f_control[row]*
                   catch_prop, size=phi))}
  ))
  post_preds_intervention_treatment= unlist(lapply(1:nrow(posterior), function(row) {
    mean(rnbinom(nrow(abundance_intervention_treatment),
                 mu = posterior$f_treatment[row]*
                   catch_prop, size=phi))}
  ))
  
  # add to vector
  ab_control_post_med[i] = median(post_preds_baseline) 
  ab_control_post_lower[i] =  quantile(post_preds_baseline, .025) 
  ab_control_post_upper[i] =  quantile(post_preds_baseline, .975) 
  
  
  ab_int_control_post_med[i] = median(post_preds_intervention_control) 
  ab_int_control_post_lower[i] =  quantile(post_preds_intervention_control, .025) 
  ab_int_control_post_upper[i] =  quantile(post_preds_intervention_control, .975) 
  
  ab_int_treatment_post_med[i] = median(post_preds_intervention_treatment) 
  ab_int_treatment_post_lower[i] =  quantile(post_preds_intervention_treatment, .025) 
  ab_int_treatment_post_upper[i] =  quantile(post_preds_intervention_treatment, .975) 
  
}

# plot vs the true values
true_ab_control = rep(NA, length(result))
for(i in 1:length(result)){
  true_ab_control[i] = mean(simulated_datasets[[i]][[3]][[1]])
}


true_ab_int_control = rep(NA, length(result))
for(i in 1:length(result)){
  true_ab_int_control[i] =  mean(simulated_datasets[[i]][[3]][[2]])
}


true_ab_int_treatment = rep(NA, length(result))
for(i in 1:length(result)){
  true_ab_int_treatment[i] =  mean(simulated_datasets[[i]][[3]][[3]])
}

# plot the true values vs the actual
setwd("~/Dropbox/DengueEntomologicalEffects/Analyses/BayesianTools/results/plots")
png("PPD_abundance_7.png", height=800, width=800)
par(mfrow=c(2,2), mar=c(5,5,4,3), cex=1.1)

# baseline
plot(0, type="n", xlim=c(0, 1500), ylim=c(0, 1500), xlab="Simulated Data", 
     ylab="Posterior Prediction", main="Baseline")
points(true_ab_control, ab_control_post_med)
abline(a=0, b=1)
segments(x0=true_ab_control, y0=ab_control_post_lower, x1=true_ab_control,
         y1=ab_control_post_upper)

# intervention control
plot(0, type="n", xlim=c(0, 1500), ylim=c(0, 1500), xlab="Simulated Data", 
     ylab="Posterior Prediction", main="Intervention Control")
points(true_ab_int_control, ab_int_control_post_med)
abline(a=0, b=1)
segments(x0=true_ab_int_control, y0=ab_int_control_post_lower, x1=true_ab_int_control,
         y1=ab_int_control_post_upper)

# intervention treatment
plot(0, type="n", xlim=c(0, 600), ylim=c(0, 600), xlab="Simulated Data", 
     ylab="Posterior Prediction", main="Intervention Treatment")
points(true_ab_int_treatment, ab_int_treatment_post_med)
abline(a=0, b=1)
segments(x0=true_ab_int_treatment, y0=ab_int_treatment_post_lower, x1=true_ab_int_treatment,
         y1=ab_int_treatment_post_upper)
dev.off()

#### PPD seroconversions ####

sc_control_post_med = rep(NA, length(result))
sc_control_post_lower = rep(NA, length(result))
sc_control_post_upper = rep(NA, length(result))

sc_treatment_post_med = rep(NA, length(result))
sc_treatment_post_lower = rep(NA, length(result))
sc_treatment_post_upper = rep(NA, length(result))

for(i in 1:length(result)){
  print(i)
  posterior = result[[i]][[1]]
  
  # use each of these to predict a new data point
  
  post_preds_baseline = unlist(lapply(1:nrow(posterior), function(row) {
    rbinom(1, size=nrow(epi[
      epi$Cluster_Allocation == "C",]),
      prob =
        1 - exp(-posterior$foi_control[row] * epi[
          epi$Cluster_Allocation == "C",]$followup_time*365))}
  ))
  post_preds_intervention = unlist(lapply(1:nrow(posterior), function(row) {
    rbinom(1, size=nrow(epi[
      epi$Cluster_Allocation == "T",]),
      prob =
        1 - exp(-posterior$foi_treatment[row] * epi[
          epi$Cluster_Allocation == "T",]$followup_time*365))}
  ))
  
  # add to vector
  sc_control_post_med[i] = median(post_preds_baseline) 
  sc_control_post_lower[i] =  quantile(post_preds_baseline, .025) 
  sc_control_post_upper[i] =  quantile(post_preds_baseline, .975) 
  
  
  sc_treatment_post_med[i] = median(post_preds_intervention) 
  sc_treatment_post_lower[i] =  quantile(post_preds_intervention, .025) 
  sc_treatment_post_upper[i] =  quantile(post_preds_intervention, .975) 
  
}

# plot vs the true values
true_sc_control = rep(NA, length(result))
for(i in 1:length(result)){
  true_sc_control[i] = sum(simulated_datasets[[i]][[1]][[1]]$outcome)
}


true_sc_treatment = rep(NA, length(result))
for(i in 1:length(result)){
  true_sc_treatment[i] = sum(simulated_datasets[[i]][[1]][[2]]$outcome)
}

# plot the true values vs the actual
setwd("~/Dropbox/DengueEntomologicalEffects/Analyses/BayesianTools/results/plots")
png("PPD_SC_7.png", height=800, width=800)
par(mfrow=c(1,2), mar=c(5,5,4,3), cex=1.1)

# baseline
plot(0, type="n", xlim=c(0,200), ylim=c(0,200), xlab="Simulated Data", 
     ylab="Posterior Prediction", main="Control")
points(true_sc_control, sc_control_post_med)
abline(a=0, b=1)
segments(x0=true_sc_control, y0=sc_control_post_lower, x1=true_sc_control,
         y1=sc_control_post_upper)

# intervention control
plot(0, type="n", xlim=c(0,200), ylim=c(0,200), xlab="Simulated Data", 
     ylab="Posterior Prediction", main="Treatment")
points(true_sc_treatment, sc_treatment_post_med)
abline(a=0, b=1)
segments(x0=true_sc_treatment, y0=sc_treatment_post_lower, x1=true_sc_treatment,
         y1=sc_treatment_post_upper)

dev.off()