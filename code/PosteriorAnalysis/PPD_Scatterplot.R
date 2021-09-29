library(tidyverse)
#### load samples, functions, trial data ####



load("~/Dropbox/DengueEntomologicalEffects/results/Posteriors/SequentialMonteCarlo.RData")

source("~/Dropbox/DengueEntomologicalEffects/code/functions.R")

epi = read.csv(file="~/Dropbox/DengueEntomologicalEffects/data/epi.csv")
parity_baseline = read.csv(file="~/Dropbox/DengueEntomologicalEffects/data/parity_baseline.csv")
parity_intervention = read.csv(file="~/Dropbox/DengueEntomologicalEffects/data/parity_intervention.csv")
load("~/Dropbox/DengueEntomologicalEffects/data/treatedclusters.RData")
abundance_baseline = read.csv(file="~/Dropbox/DengueEntomologicalEffects/data/abundance_baseline.csv")
abundance_intervention_control = read.csv(file="~/Dropbox/DengueEntomologicalEffects/data/abundance_intervention_control.csv")
abundance_intervention_treatment = read.csv(file="~/Dropbox/DengueEntomologicalEffects/data/abundance_intervention_treatment.csv")
bloodmeal_baseline = read.csv(file="~/Dropbox/DengueEntomologicalEffects/data/bloodmeal_baseline.csv")
bloodmeal_intervention = read.csv(file="~/Dropbox/DengueEntomologicalEffects/data/bloodmeal_intervention.csv")

untreated.clusters = seq(1,26)[-treated.clusters]

C = 0.45
p = 5.2

samples=head(samples, 1000)

#### modified biting and death rates ####

samples$gc_control = rep(NA, nrow(samples))
samples$gc_treatment = rep(NA, nrow(samples))

mosquito.locations.control = matrix(rep(NA, nrow(samples)*9), nrow=nrow(samples),
                                    ncol=9)

for(i in 1:nrow(samples)){
  mosquito.locations.control[i,] = with(samples[i,],  location_probability(af, ap, qu, tau, d, rho,
                                                                           qt, alpha, C=0))
} 

mosquito.locations.treatment = matrix(rep(NA, nrow(samples)*9), nrow=nrow(samples),
                                      ncol=9)

for(i in 1:nrow(samples)){
  mosquito.locations.treatment[i,] = with(samples[i,],  location_probability(af, ap,
                                                                             qu, tau,
                                                                             d, rho,
                                                                             qt, alpha,
                                                                             C=.47))
} 

samples$pi.U.control =  mosquito.locations.control[,1] +
  mosquito.locations.control[,2]+mosquito.locations.control[,3]
samples$pi.tau.control =  mosquito.locations.control[,4] +
  mosquito.locations.control[,5]+mosquito.locations.control[,6]
samples$pi.T.control =  mosquito.locations.control[,7] +
  mosquito.locations.control[,8]+mosquito.locations.control[,9]

samples$pi.U.treated =  mosquito.locations.treatment[,1] +
  mosquito.locations.treatment[,2]+mosquito.locations.treatment[,3]
samples$pi.tau.treated =  mosquito.locations.treatment[,4] +
  mosquito.locations.treatment[,5]+mosquito.locations.treatment[,6]
samples$pi.T.treated =  mosquito.locations.treatment[,7] +
  mosquito.locations.treatment[,8]+mosquito.locations.treatment[,9]

for(i in 1:nrow(samples)){
  samples$gc_control[i] = with(samples[i,], gc_calc(mu, qu, qt, gu, gt, gtau, tau,
                                                    rho, C=0, 
                                                    pi.tau=pi.tau.control, pi.U.control,
                                                    pi.T.control))
  samples$gc_treatment[i] = with(samples[i,], gc_calc(mu, qu, qt, gu, gt, gtau, tau,
                                                      rho, C=.45,
                                                      pi.tau=pi.tau.treated, pi.U.treated,
                                                      pi.T.treated))
}

samples$ac_control = samples$ac_treatment = rep(NA, nrow(samples))
for(i in 1:nrow(samples)){
  samples$ac_control[i] = with(samples[i,], ac_calc(af, ap, rho, C=0, alpha, qt, qu,
                                                    tau))
  samples$ac_treatment[i] = with(samples[i,], ac_calc(af, ap, rho, C=.45, alpha, qt, qu,
                                                      tau))
}


#### calculate incidence, parity, mosquito abundance, and blood meal status from posterior ####

parity.c = rep(NA, nrow(samples))
parity.t = rep(NA, nrow(samples))

for(i in 1:nrow(samples)){
  parity.c[i] = parity(samples[i, "ac_control"],samples[i, "gc_control"])
  parity.t[i] = parity(samples[i, "ac_treatment"], samples[i, "gc_treatment"])
}

foi.matrix.c = foi.matrix.t = 
  matrix(rep(NA, 26*nrow(samples)), ncol=26)

foi.matrix.c = foi.matrix.t = 
  matrix(rep(NA, 26*nrow(samples)), ncol=26)

for(i in 1:nrow(samples)){
  print(i)
  for(j in 1:26){
    foi.matrix.c[i,j] = with(samples[i,], force_of_infection(b, lambda, ac_control, 
                                                             c, get(paste0("X_", j)),
                                                             gc_control, n))
    foi.matrix.t[i,j] = with(samples[i,], force_of_infection(b, lambda, ac_treatment, 
                                                             c, get(paste0("X_", j)),
                                                             gc_treatment, n))
  }
}

prob.untreated.c = prob.untreated.t = 
  prob.treated.t = rep(NA,nrow(samples))


# modeled mosquito to human ratio, control and treatment
m_control = samples$lambda / samples$gc_control
m_treatment = samples$lambda / samples$gc_treatment


prob.unfed.c = prob.partfed.c=prob.fullfed.c=
  prob.unfed.t=
  prob.partfed.t=
  prob.fullfed.t = prob.untreated.c = prob.untreated.t = 
  prob.treated.t = rep(NA,nrow(samples))


prob.fullyfed.c = mosquito.locations.control[,1] / samples$pi.U.control
prob.partfed.c =  mosquito.locations.control[,2] / samples$pi.U.control
prob.unfed.c =  mosquito.locations.control[,3] / samples$pi.U.control

prob.untreated.c = mosquito.locations.control[,1] + mosquito.locations.control[,2]+
  mosquito.locations.control[,3]

prob.fullyfed.t = mosquito.locations.treatment[,7] / samples$pi.T.treated
prob.partfed.t =  mosquito.locations.treatment[,8] / samples$pi.T.treated
prob.unfed.t =  mosquito.locations.treatment[,9] / samples$pi.T.treated

prob.treated.t = mosquito.locations.treatment[,7] + mosquito.locations.treatment[,8]+
  mosquito.locations.treatment[,9]

expected_mosquitoes_control = p * m_control * prob.untreated.c * samples$catch_prop
expected_mosquitoes_treatment = p * m_treatment * prob.treated.t * samples$catch_prop


#### Seroconversion Data ####

png("~/Dropbox/DengueEntomologicalEffects/results/ManuscriptPlots/PPDAllScatter.png",
    width=800, height=800)

par(mfrow=c(2,2), cex=1.2)

plot(0, type="n", xlim=c(0, 1), ylim=c(0, 1),
     xlab="Predicted proportion seroconverting",
     ylab="Proportion Seroconverting")

for(i in 1:length(untreated.clusters)){
  
  # extract seroconversion data
  sc = epi %>%
    filter(Cluster == untreated.clusters[i]) %>%
    pull(outcome)
  # total number of seroconversions
  true.sc = sum(sc)
  
  # followup time
  fut = epi %>%
    filter(Cluster == untreated.clusters[i]) %>%
    pull(followup_time)
  
  # number of individuals
  n.ind = length(sc)
  
  
  # find the probability that each individual did not
  # seroconvert, given the posterior draw
  prob_didnt_seroconvert = matrix(rep(NA, nrow(samples)*n.ind),
                                  nrow=nrow(samples), ncol=n.ind)
  for(j in 1:nrow(samples)){
    for(k in 1:n.ind){
      prob_didnt_seroconvert[j,k] = exp(-foi.matrix.c[j,untreated.clusters[i]]*fut[k])
    }
  }
  
  # Randomly draw values from these for if they didnt seroconvert
  didnt_seroconvert_sim = matrix(rep(NA, nrow(samples)*n.ind),
                                 nrow=nrow(samples), ncol=n.ind)
  for(j in 1:nrow(samples)){
    for(k in 1:n.ind){
      didnt_seroconvert_sim[j,k] = rbinom(1,
                                          size=1,
                                          prob_didnt_seroconvert[j,k])
    }
  }
  
  # change to if they did seroconvert
  did_seroconvert_sim = 1-didnt_seroconvert_sim
  
  # find total number seroconverting in the cluster by
  # summing over rows
  points(mean(rowSums(did_seroconvert_sim)/n.ind), sum(sc)/n.ind, pch=8, col=rgb(1,0,0,alpha=.5))
  
  
  
}
abline(a=0, b=1)

for(i in 1:length(treated.clusters)){
  
  # extract seroconversion data
  sc = epi %>%
    filter(Cluster == treated.clusters[i]) %>%
    pull(outcome)
  # total number of seroconversions
  true.sc = sum(sc)
  
  # followup time
  fut = epi %>%
    filter(Cluster == treated.clusters[i]) %>%
    pull(followup_time)
  
  # number of individuals
  n.ind = length(sc)
  
  
  # find the probability that each individual did not
  # seroconvert, given the posterior draw
  prob_didnt_seroconvert = matrix(rep(NA, nrow(samples)*n.ind),
                                  nrow=nrow(samples), ncol=n.ind)
  for(j in 1:nrow(samples)){
    for(k in 1:n.ind){
      prob_didnt_seroconvert[j,k] = exp(-foi.matrix.t[j,untreated.clusters[i]]*fut[k])
    }
  }
  
  # Randomly draw values from these for if they didnt seroconvert
  didnt_seroconvert_sim = matrix(rep(NA, nrow(samples)*n.ind),
                                 nrow=nrow(samples), ncol=n.ind)
  for(j in 1:nrow(samples)){
    for(k in 1:n.ind){
      didnt_seroconvert_sim[j,k] = rbinom(1,
                                          size=1,
                                          prob_didnt_seroconvert[j,k])
    }
  }
  
  # change to if they did seroconvert
  did_seroconvert_sim = 1-didnt_seroconvert_sim
  
  # find total number seroconverting in the cluster by
  # summing over rows
  points(mean(rowSums(did_seroconvert_sim)/n.ind), sum(sc)/n.ind, pch=19, col=rgb(0,1,0,alpha=.5))
  
  
  
}


#### Parity ####

plot(0, type="n", xlim=c(0, 1), ylim=c(0, 1),
     xlab="Predicted proportion parous",
     ylab="True proportion parous")
     
     
for(i in 1:length(treated.clusters)){
  cluster = treated.clusters[i]
  
  sample_size = parity_baseline %>%
    filter(Cluster  == cluster) %>%
    pull(total_caught)
  number_parous = parity_baseline %>%
    filter(Cluster  == cluster) %>%
    pull(parous)
  
  points(jitter(mean(parity.c)), number_parous/sample_size, pch=19, col=rgb(1,0,0,alpha=.5))
  
}



for(i in 1:length(untreated.clusters)){
  cluster = untreated.clusters[i]
  
  sample_size = parity_baseline %>%
    filter(Cluster  == cluster) %>%
    pull(total_caught)
  number_parous = parity_baseline %>%
    filter(Cluster  == cluster) %>%
    pull(parous)
  
  points(jitter(mean(parity.c)), number_parous/sample_size, pch=8, col=rgb(1,0,0,alpha=.5))
  
}

for(i in 1:length(treated.clusters)){
  cluster = treated.clusters[i]
  
  sample_size = parity_intervention %>%
    filter(Cluster  == cluster) %>%
    pull(total_caught)
  number_parous = parity_intervention %>%
    filter(Cluster  == cluster) %>%
    pull(parous)
  
  points(jitter(mean(parity.t)), number_parous/sample_size, pch=19, col=rgb(0,1,0,alpha=.5))
  
}
abline(0,1)



for(i in 1:length(untreated.clusters)){
  cluster = untreated.clusters[i]
  
  sample_size = parity_intervention %>%
    filter(Cluster  == cluster) %>%
    pull(total_caught)
  number_parous = parity_intervention %>%
    filter(Cluster  == cluster) %>%
    pull(parous)
  
  points(jitter(mean(parity.c)), number_parous/sample_size, pch=8, col=rgb(0,1,0,alpha=.5))
  
}

legend("topleft", c("Control Clusters", "Treated Clusters", "Baseline Period", "Intervention Period"),
       col = c("black", "black", rgb(1,0,0,alpha=.5), rgb(0,1,0,alpha=.5)),
       pch=c(8, 19, 15, 15))

#### mosquito aspiration ####

plot(0, type="n", xlim=c(0.2, 20), ylim=c(0.2, 0.8),
     xlab="Predicted catch number",
     ylab="Actual catch number")

for(i in 1:length(untreated.clusters)){
  
  abundance_data = abundance_baseline %>%
    filter(Cluster == untreated.clusters[i]) %>%
    pull(total_caught)
  
  sim_asp_means = rep(NA, nrow(samples))
  
  for(j in 1:nrow(samples)){
    sim_asp_means[j]= mean(rnbinom(length(abundance_data),
                                   mu=expected_mosquitoes_control[j],
                                   size=samples$phi[j]))
  }
  
  points(jitter(mean(sim_asp_means)), jitter(mean(abundance_data)), pch=8, col=rgb(1,0,0,alpha=.5))
  
}

for(i in 1:length(treated.clusters)){
  
  abundance_data = abundance_baseline %>%
    filter(Cluster == treated.clusters[i]) %>%
    pull(total_caught)
  
  sim_asp_means = rep(NA, nrow(samples))
  
  for(j in 1:nrow(samples)){
    sim_asp_means[j]= mean(rnbinom(length(abundance_data),
                                   mu=expected_mosquitoes_control[j],
                                   size=samples$phi[j]))
  }
  
  points(jitter(mean(sim_asp_means)), jitter(mean(abundance_data)), pch=19, col=rgb(1,0,0,alpha=.5))
  
}

for(i in 1:length(untreated.clusters)){
  
  abundance_data = abundance_intervention_control %>%
    filter(Cluster == untreated.clusters[i]) %>%
    pull(total_caught)
  
  sim_asp_means = rep(NA, nrow(samples))
  
  for(j in 1:nrow(samples)){
    sim_asp_means[j]= mean(rnbinom(length(abundance_data),
                                   mu=expected_mosquitoes_control[j],
                                   size=samples$phi[j]))
  }
  
  points(jitter(mean(sim_asp_means)), jitter(mean(abundance_data)), pch=8, col=rgb(0,1,0,alpha=.5))
  
}

for(i in 1:length(treated.clusters)){
  print(i)
  
  abundance_data = abundance_intervention_treatment %>%
    filter(Cluster == treated.clusters[i]) %>%
    pull(total_caught)
  
  sim_asp_means = rep(NA, nrow(samples))
  
  for(j in 1:nrow(samples)){
    sim_asp_means[j]= mean(rnbinom(length(abundance_data),
                                   mu=expected_mosquitoes_treatment[j],
                                   size=samples$phi[j]))
  }
  
  points(jitter(mean(sim_asp_means)), jitter(mean(abundance_data)), pch=19, col=rgb(0,1,0,alpha=.5))
  
}

abline(a=0, b=1)


#### Bloodmeal status ####

plot(0, type="n", xlim=c(0.2, 0.7), ylim=c(0.2, 0.7),
     xlab="Predicted proportion unfed",
     ylab="Actual proportion unfed")

for(i in 1:length(untreated.clusters)){
  
  unfed_prop = bloodmeal_baseline %>%
    filter(Cluster == untreated.clusters[i]) %>%
    mutate(prop_unfed = total_empty / total_evaluated) %>%
    pull(prop_unfed)
  
  prob.unfed.c
  
  points(jitter(mean(prob.unfed.c)), unfed_prop, pch=8, col=rgb(1,0,0,alpha=.5))
  
}

for(i in 1:length(treated.clusters)){
  
  unfed_prop = bloodmeal_baseline %>%
    filter(Cluster == treated.clusters[i]) %>%
    mutate(prop_unfed = total_empty / total_evaluated) %>%
    pull(prop_unfed)
  
  prob.unfed.c
  
  points(jitter(mean(prob.unfed.c)), unfed_prop, pch=19, col=rgb(1,0,0,alpha=.5))
  
}

for(i in 1:length(untreated.clusters)){
  
  unfed_prop = bloodmeal_intervention %>%
    filter(Cluster == untreated.clusters[i]) %>%
    mutate(prop_unfed = total_empty / total_evaluated) %>%
    pull(prop_unfed)
  
  prob.unfed.c
  
  points(jitter(mean(prob.unfed.c)), unfed_prop, pch=8, col=rgb(0,1,0,alpha=.5))
  
}

for(i in 1:length(treated.clusters)){
  
  unfed_prop = bloodmeal_intervention %>%
    filter(Cluster == treated.clusters[i]) %>%
    mutate(prop_unfed = total_empty / total_evaluated) %>%
    pull(prop_unfed)
  
  points(jitter(mean(prob.unfed.t)), unfed_prop, pch=19, col=rgb(0,1,0,alpha=.5))
  
}

abline(a=0, b=1)

dev.off()

#### Comparison of proportion unfed mosquitoes and parity rate ####
par(mfrow=c(1,1), mar=c(5,5,2,2))
plot(0, type="n", xlab="Proportion mosquitoes unfed",
     ylab="Proportion parous", xlim=c(0,1), ylim=c(0,1))

for(i in 1:length(untreated.clusters)){
  print(i)
  unfed_prop = bloodmeal_baseline %>%
    filter(Cluster == untreated.clusters[i]) %>%
    mutate(prop_unfed = total_empty / total_evaluated) %>%
    pull(prop_unfed)
  
  prop_parous = parity_baseline %>%
    filter(Cluster  == untreated.clusters[i]) %>%
    mutate(prop_parous = parous/total_caught) %>%
    pull(prop_parous)
  
  points(unfed_prop, prop_parous, pch=19, col="red")
}

for(i in 1:length(treated.clusters)){
  print(i)
  unfed_prop = bloodmeal_baseline %>%
    filter(Cluster == treated.clusters[i]) %>%
    mutate(prop_unfed = total_empty / total_evaluated) %>%
    pull(prop_unfed)
  
  prop_parous = parity_baseline %>%
    filter(Cluster  == treated.clusters[i]) %>%
    mutate(prop_parous = parous/total_caught) %>%
    pull(prop_parous)
  
  points(unfed_prop, prop_parous, pch=19, col="blue")
}

for(i in 1:length(untreated.clusters)){
  print(i)
  unfed_prop = bloodmeal_intervention %>%
    filter(Cluster == untreated.clusters[i]) %>%
    mutate(prop_unfed = total_empty / total_evaluated) %>%
    pull(prop_unfed)
  
  prop_parous = parity_intervention %>%
    filter(Cluster  == untreated.clusters[i]) %>%
    mutate(prop_parous = parous/total_caught) %>%
    pull(prop_parous)
  
  points(unfed_prop, prop_parous, pch=19, col="green")
}

for(i in 1:length(treated.clusters)){
  print(i)
  unfed_prop = bloodmeal_intervention %>%
    filter(Cluster == treated.clusters[i]) %>%
    mutate(prop_unfed = total_empty / total_evaluated) %>%
    pull(prop_unfed)
  
  prop_parous = parity_intervention %>%
    filter(Cluster  == treated.clusters[i]) %>%
    mutate(prop_parous = parous/total_caught) %>%
    pull(prop_parous)
  
  points(unfed_prop, prop_parous, pch=19, col="yellow")
}