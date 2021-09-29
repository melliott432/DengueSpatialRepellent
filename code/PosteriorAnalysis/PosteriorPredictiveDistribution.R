library(tidyverse)
library(MCMCpack)
#### load samples, functions, trial data ####



load("~/Dropbox/DengueEntomologicalEffects/results/Posteriors/SMC_Taper4.RData")

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
  
prob.fullyfed.t = mosquito.locations.treatment[,7] / samples$pi.T.treated
prob.partfed.t =  mosquito.locations.treatment[,8] / samples$pi.T.treated
prob.unfed.t =  mosquito.locations.treatment[,9] / samples$pi.T.treated



#### seroconversion PPD ####

png("~/Dropbox/DengueEntomologicalEffects/results/ManuscriptPlots/PPDAll.png",
    height=1000, width=1000)
par(mfrow=c(1,1), cex=1.2, oma=c(1, 2, 3, 0), mar=c(3, 3, 3, 1))

plot(0, type = "n", xlim=c(1,26), ylim=c(-.1,1), 
     main="Seroconversions", xaxt="n", xlab="", 
     xaxs="i", yaxs="i", xaxt="n",
     ylab="Proportion Seroconverting")
axis(1, at=8,
     labels="Control Clusters",
     tick=F)

axis(1, at=20, labels="Treated Clusters", tick=F)
rect.colors = c(rep(c("lightgrey", "white"), 14))
for(i in 1:28){
  rect(i, -2, i+1, 600, col=rect.colors[i],
       border=NA)
}
abline(v=14)

# posterior for bernoulli draw
beta_bernoulli = function(x){
  if(x == 0){
    posterior_draw = rbinom(1000, size=1, rbeta(1, 1, 2))
  } else{
    posterior_draw = rbinom(1000, size=1, rbeta(1, 2, 1))
  }
  return(posterior_draw)
}


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
  boxplot(rowSums(did_seroconvert_sim)/n.ind, at=i+.25, add=T)
  
  # add data
  points(i+.75, sum(sc)/n.ind, col="red", pch=19)
  
  
}

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
      prob_didnt_seroconvert[j,k] = exp(-foi.matrix.t[j,treated.clusters[i]]*fut[k])
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
  boxplot(rowSums(did_seroconvert_sim)/ n.ind, at=i+.25+13, add=T)
  
  # add a point for data
  points(i+.75+13, sum(sc)/n.ind, col="red", pch=19)
  
}

legend("topright", c("Model \n Prediction", "Trial Data"),
       col = c("black", "red"),
       pch=19)

#### parity baseline ####

plot(0, type = "n", xlim=c(2,26), ylim=c(0,1), 
     main="Parity, Baseline", xaxt="n", xlab="", 
     xaxs="i", yaxs="i", xaxt="n", ylab="Proportion Parous")
axis(1, at=8,
     labels="Control Clusters",
     tick=F)

axis(1, at=20, labels="Treated Clusters", tick=F)
rect.colors = c(rep(c("lightgrey", "white"), 14))
for(i in 1:28){
  rect(i, -.5, i+1, 600, col=rect.colors[i],
       border=NA)
}
abline(v=14)

for(i in 1:length(treated.clusters)){
  cluster = treated.clusters[i]
  
  
  sample_size = parity_baseline %>%
    filter(Cluster  == i) %>%
    pull(total_caught)
  number_parous = parity_baseline %>%
    filter(Cluster  == i) %>%
    pull(parous)
  
  # model data
  boxplot(rbinom(1000, size=sample_size, 
                 prob = parity.c)/sample_size, add=T, 
          at=i+.25+13, col="blue")
  
  # data from trial
  data_draws = rbeta(1000,
                     1+number_parous, 
                     1 + sample_size - number_parous)
  boxplot(rbinom(1000, size=sample_size, 
                 prob = data_draws)/sample_size, add=T, 
          at=i+13+.75, col="red")
  
}



for(i in 1:length(seq(1,26)[-treated.clusters])){
  cluster = seq(1,26)[-treated.clusters][i]
  
  
  sample_size = parity_baseline %>%
    filter(Cluster  == i) %>%
    pull(total_caught)
  number_parous = parity_baseline %>%
    filter(Cluster  == i) %>%
    pull(parous)
  
  # model data
  boxplot(rbinom(1000, size=sample_size, 
                 prob = parity.c)/sample_size, add=T, 
          at=i+.25, col="blue")
  
  # data from trial
  data_draws = rbeta(1000,
                     1+number_parous, 
                     1 + sample_size - number_parous)
  boxplot(rbinom(1000, size=sample_size, 
                 prob = data_draws)/sample_size, add=T, 
          at=i+.75, col="red")
  
}
legend("topleft", c("Model", "Data"), col=c("blue", "red"), pch=19)


#### Parity Intervention ####

plot(0, type = "n", xlim=c(2,26), ylim=c(0,1), 
     main="Parity, Intervention", xaxt="n", xlab="", 
     xaxs="i", yaxs="i", xaxt="n", ylab="")
rect.colors = c(rep(c("lightgrey", "white"), 14))
for(i in 1:28){
  rect(i, -.5, i+1, 600, col=rect.colors[i],
       border=NA)
}
axis(1, at=8,
     labels="Control Clusters",
     tick=F)

axis(1, at=20, labels="Treated Clusters", tick=F)
abline(v=14)


for(i in 1:length(treated.clusters)){
  cluster = treated.clusters[i]
  
  
  sample_size = parity_intervention %>%
    filter(Cluster  == i) %>%
    pull(total_caught)
  number_parous = parity_intervention %>%
    filter(Cluster  == i) %>%
    pull(parous)
  
  # model data
  boxplot(rbinom(1000, size=sample_size, 
                 prob = parity.t)/sample_size, add=T, 
          at=i+.25+13, col="blue")
  
  # data from trial
  data_draws = rbeta(1000,
                     1+number_parous, 
                     1 + sample_size - number_parous)
  boxplot(rbinom(1000, size=sample_size, 
                 prob = data_draws)/sample_size, add=T, 
          at=i+13+.75, col="red")
  
}



for(i in 1:length(seq(1,26)[-treated.clusters])){
  cluster = seq(1,26)[-treated.clusters][i]
  
  
  sample_size = parity_intervention %>%
    filter(Cluster  == i) %>%
    pull(total_caught)
  number_parous = parity_intervention %>%
    filter(Cluster  == i) %>%
    pull(parous)
  
  # model data
  boxplot(rbinom(1000, size=sample_size, 
                 prob = parity.c)/sample_size, add=T, 
          at=i+.25, col="blue")
  
  # data from trial
  data_draws = rbeta(1000,
                     1+number_parous, 
                     1 + sample_size - number_parous)
  boxplot(rbinom(1000, size=sample_size, 
                 prob = data_draws)/sample_size, add=T, 
          at=i+.75, col="red")
  
}

#### Mosquito Aspirations ####



# expected number of mosquitoes present in aspiration
expected_mosquitoes_control = p * m_control * samples$pi.U.control
expected_mosquitoes_treatment =  p * m_treatment * samples$pi.T.treated 

plot(0, type = "n", xlim=c(1,26), ylim=c(-.1,1.5), 
     main="Mosquito Abundance, \nBaseline", xaxt="n", xlab="", 
     xaxs="i", yaxs="i", xaxt="n", ylab="Mean Mosquitoes per Aspiration")
axis(1, at=8,
     labels="Control Clusters",
     tick=F)

axis(1, at=20, labels="Treated Clusters", tick=F)
rect.colors = c(rep(c("lightgrey", "white"), 14))
for(i in 1:28){
  rect(i, -.5, i+1, 600, col=rect.colors[i],
       border=NA)
}
abline(v=14)

for(i in 1:length(untreated.clusters)){
  
  abundance_data = abundance_baseline %>%
    filter(Cluster == untreated.clusters[i]) %>%
    pull(total_caught)
  points(i+.75, mean(abundance_data), col="red", pch=19)
  
  sim_asp_means = rep(NA, nrow(samples))
  
  for(j in 1:nrow(samples)){
    sim_asp_means[j]= mean(rnbinom(length(abundance_data),
                           mu=expected_mosquitoes_control[j]*samples$catch_prop,
                           size=samples$phi[j]))
  }
  
  boxplot(sim_asp_means, add=T, at = i+.25)

}

for(i in 1:length(treated.clusters)){
  
  abundance_data = abundance_baseline %>%
    filter(Cluster == untreated.clusters[i]) %>%
    pull(total_caught)
  points(i+13+.75, mean(abundance_data), col="red", pch=19)
  
  sim_asp_means = rep(NA, nrow(samples))
  
  for(j in 1:nrow(samples)){
    sim_asp_means[j]= mean(rnbinom(length(abundance_data),
                                   mu=expected_mosquitoes_control[j]*samples$catch_prop,
                                   size=samples$phi[j]))
  }
  
  boxplot(sim_asp_means, add=T, at = i+13+.25)
}


# Abundance Intervention

plot(0, type = "n", xlim=c(1,26), ylim=c(-.1,1.5), 
     main="Mosquito Abundance, \nIntervention", xaxt="n", xlab="", 
     xaxs="i", yaxs="i", xaxt="n", ylab="Mean Mosquitoes per Aspiration")
rect.colors = c(rep(c("lightgrey", "white"), 14))
for(i in 1:28){
  rect(i, -.5, i+1, 600, col=rect.colors[i],
       border=NA)
}
axis(1, at=8,
     labels="Control Clusters",
     tick=F)

axis(1, at=20, labels="Treated Clusters", tick=F)
abline(v=14)

for(i in 1:length(untreated.clusters)){
  
  abundance_data = abundance_intervention_control %>%
    filter(Cluster == untreated.clusters[i]) %>%
    pull(total_caught)
  points(i+.75, mean(abundance_data), col="red", pch=19)
  
  sim_asp_means = rep(NA, nrow(samples))
  
  for(j in 1:nrow(samples)){
    sim_asp_means[j]= mean(rnbinom(length(abundance_data),
                                   mu=expected_mosquitoes_control[j]*samples$catch_prop,
                                   size=samples$phi[j]))
  }
  
  boxplot(sim_asp_means, add=T, at = i+.25)

  
}

for(i in 1:length(treated.clusters)){
  
  abundance_data = abundance_intervention_treatment %>%
    filter(Cluster == treated.clusters[i]) %>%
    pull(total_caught)
  points(i+13+.75, mean(abundance_data), col="red", pch=19)
  
  sim_asp_means = rep(NA, nrow(samples))
  
  for(j in 1:nrow(samples)){
    sim_asp_means[j]= mean(rnbinom(length(abundance_data),
                                   mu=expected_mosquitoes_treatment[j]*samples$catch_prop,
                                   size=samples$phi[j]))
  }
  
  boxplot(sim_asp_means, add=T, at = i+13+.25)
  
}


##### Blood feeding baseline ####

plot(0, type="n", xlim=c(1, 26), ylim=c(0.1,0.9), 
     main="Bloodmeal, Baseline", xlab="", xaxt="n")
axis(side =1, 1:13, labels = paste0("Cluster", 1:13))

axis(1, at=8,
     labels="Control Clusters",
     tick=F)

axis(1, at=20, labels="Treated Clusters", tick=F)

for(i in 1:length(untreated.clusters)){
  unfed = bloodmeal_baseline %>%
    filter(Cluster == untreated.clusters[i]) %>%
    mutate(propn_unfed = total_empty / total_evaluated) %>%
    pull(propn_unfed)
  
  partial = bloodmeal_baseline %>%
    filter(Cluster == untreated.clusters[i]) %>%
    mutate(propn_partfed = total_half / total_evaluated) %>%
    pull(propn_partfed)
  
  fully = bloodmeal_baseline %>%
    filter(Cluster == untreated.clusters[i]) %>%
    mutate(propn_fullfed = total_full / total_evaluated) %>%
    pull(propn_fullfed)
  
  total_evaluated = bloodmeal_baseline %>%
    filter(Cluster == untreated.clusters[i]) %>%
   pull(total_evaluated)
  
  un = total_evaluated*unfed
  part = total_evaluated*partial
  full = total_evaluated*fully
  
  polygon(x=c(i, i, i-.25, i-.25), c(0, unfed, unfed, 0), col="red")
  polygon(x=c(i, i, i-.25, i-.25), c(unfed, partial+unfed,partial+unfed, unfed),
          col="blue")
  polygon(x=c(i, i, i-.25, i-.25), c(partial+unfed, partial+fully+unfed,
                                         partial+fully+unfed,
                                         partial+unfed),
          col="green")
  
  # uncertainty
  alpha.vec = c(1+un, 1+part, 1+full)
  uncertainty = rdirichlet(1000, alpha.vec)
  
  lines(x=c(i-.125, i-.125), y=c(quantile(uncertainty[,1], 0.025),
                               quantile(uncertainty[,1], 0.975)))
  lines(x=c(i-.125, i-.125), y=c(quantile(uncertainty[,2], 0.025)+unfed,
                               quantile(uncertainty[,2], 0.975)+unfed))
  lines(x=c(i-.125, i-.125), y=c(quantile(uncertainty[,3], 0.025)+
                                 unfed+partial,
                               quantile(uncertainty[,3], 0.975)+
                                 unfed+partial))
  
  # plot the mode-predicted proportions
  mean.unfed = mean(prob.unfed.c)
  mean.partial = mean(prob.partfed.c)
  mean.fully = mean(prob.fullyfed.c)
  polygon(x=c(i+.25, i+.25, i, i), y=c(0, mean.unfed, mean.unfed, 0), 
          col=alpha("red", .5))
  polygon(x=c(i+.25, i+.25, i, i), c(mean.unfed, mean.partial+
                                     mean.unfed,mean.partial+mean.unfed,
                                   mean.unfed),
          col=alpha("blue", .5))
  polygon(x=c(i+.25, i+.25, i, i), c(mean.partial+mean.unfed,
                                   mean.partial+mean.fully+mean.unfed,
                                   mean.partial+mean.fully+mean.unfed,
                                   mean.partial+mean.unfed),
          col=alpha("green", .5))
  
  
}


for(i in 1:length(treated.clusters)){
  unfed = bloodmeal_baseline %>%
    filter(Cluster == treated.clusters[i]) %>%
    mutate(propn_unfed = total_empty / total_evaluated) %>%
    pull(propn_unfed)
  
  partial = bloodmeal_baseline %>%
    filter(Cluster == treated.clusters[i]) %>%
    mutate(propn_partfed = total_half / total_evaluated) %>%
    pull(propn_partfed)
  
  fully = bloodmeal_baseline %>%
    filter(Cluster == treated.clusters[i]) %>%
    mutate(propn_fullfed = total_full / total_evaluated) %>%
    pull(propn_fullfed)
  
  total_evaluated = bloodmeal_baseline %>%
    filter(Cluster == treated.clusters[i]) %>%
    pull(total_evaluated)
  
  un = total_evaluated*unfed
  part = total_evaluated*partial
  full = total_evaluated*fully
  
  polygon(x=c(i+13, i+13, i+13-.25, i+13-.25), c(0, unfed, unfed, 0), col="red")
  polygon(x=c(i+13, i+13, i+13-.25, i+13-.25), c(unfed, partial+unfed,partial+unfed, unfed),
          col="blue")
  polygon(x=c(i+13, i+13, i+13-.25, i+13-.25), c(partial+unfed, partial+fully+unfed,
                                     partial+fully+unfed,
                                     partial+unfed),
          col="green")
  
  # uncertainty
  alpha.vec = c(1+un, 1+part, 1+full)
  uncertainty = rdirichlet(1000, alpha.vec)
  
  lines(x=c(i+13-.125, i+13-.125), y=c(quantile(uncertainty[,1], 0.025),
                               quantile(uncertainty[,1], 0.975)))
  lines(x=c(i+13-.125, i+13-.125), y=c(quantile(uncertainty[,2], 0.025)+unfed,
                               quantile(uncertainty[,2], 0.975)+unfed))
  lines(x=c(i+13-.125, i+13-.125), y=c(quantile(uncertainty[,3], 0.025)+
                                 unfed+partial,
                               quantile(uncertainty[,3], 0.975)+
                                 unfed+partial))
  
  # plot the mode-predicted proportions
  mean.unfed = mean(prob.unfed.c)
  mean.partial = mean(prob.partfed.c)
  mean.fully = mean(prob.fullyfed.c)
  polygon(x=c(i+13+.25, i+13+.25, i+13, i+13), y=c(0, mean.unfed, mean.unfed, 0), 
          col=alpha("red", .5))
  polygon(x=c(i+13+.25, i+13+.25, i+13, i+13), c(mean.unfed, mean.partial+
                                       mean.unfed,mean.partial+mean.unfed,
                                     mean.unfed),
          col=alpha("blue", .5))
  polygon(x=c(i+13+.25, i+13+.25, i+13, i+13), c(mean.partial+mean.unfed,
                                     mean.partial+mean.fully+mean.unfed,
                                     mean.partial+mean.fully+mean.unfed,
                                     mean.partial+mean.unfed),
          col=alpha("green", .5))
  
  
}
abline(v=13.5, lwd=4)

##### Mosquito feeding intervention ####
plot(0, type="n", xlim=c(1, 26), ylim=c(0.1,0.9), 
     main="Bloodmeal, Intervention", xlab="", xaxt="n")
axis(side =1, 1:13, labels = paste0("Cluster", 1:13))
axis(1, at=8,
     labels="Control Clusters",
     tick=F)

axis(1, at=20, labels="Treated Clusters", tick=F)



for(i in 1:length(untreated.clusters)){
  unfed = bloodmeal_intervention %>%
    filter(Cluster == untreated.clusters[i]) %>%
    mutate(propn_unfed = total_empty / total_evaluated) %>%
    pull(propn_unfed)
  
  partial = bloodmeal_intervention %>%
    filter(Cluster == untreated.clusters[i]) %>%
    mutate(propn_partfed = total_half / total_evaluated) %>%
    pull(propn_partfed)
  
  fully = bloodmeal_intervention %>%
    filter(Cluster == untreated.clusters[i]) %>%
    mutate(propn_fullfed = total_full / total_evaluated) %>%
    pull(propn_fullfed)
  
  total_evaluated = bloodmeal_intervention %>%
    filter(Cluster == untreated.clusters[i]) %>%
    pull(total_evaluated)
  
  un = total_evaluated*unfed
  part = total_evaluated*partial
  full = total_evaluated*fully
  
  polygon(x=c(i, i, i-.25, i-.25), c(0, unfed, unfed, 0), col="red")
  polygon(x=c(i, i, i-.25, i-.25), c(unfed, partial+unfed,partial+unfed, unfed),
          col="blue")
  polygon(x=c(i, i, i-.25, i-.25), c(partial+unfed, partial+fully+unfed,
                                     partial+fully+unfed,
                                     partial+unfed),
          col="green")
  
  # uncertainty
  alpha.vec = c(1+un, 1+part, 1+full)
  uncertainty = rdirichlet(1000, alpha.vec)
  
  lines(x=c(i-.125, i-.125), y=c(quantile(uncertainty[,1], 0.025),
                                 quantile(uncertainty[,1], 0.975)))
  lines(x=c(i-.125, i-.125), y=c(quantile(uncertainty[,2], 0.025)+unfed,
                                 quantile(uncertainty[,2], 0.975)+unfed))
  lines(x=c(i-.125, i-.125), y=c(quantile(uncertainty[,3], 0.025)+
                                   unfed+partial,
                                 quantile(uncertainty[,3], 0.975)+
                                   unfed+partial))
  
  # plot the mode-predicted proportions
  mean.unfed = mean(prob.unfed.c)
  mean.partial = mean(prob.partfed.c)
  mean.fully = mean(prob.fullyfed.c)
  polygon(x=c(i+.25, i+.25, i, i), y=c(0, mean.unfed, mean.unfed, 0), 
          col=alpha("red", .5))
  polygon(x=c(i+.25, i+.25, i, i), c(mean.unfed, mean.partial+
                                       mean.unfed,mean.partial+mean.unfed,
                                     mean.unfed),
          col=alpha("blue", .5))
  polygon(x=c(i+.25, i+.25, i, i), c(mean.partial+mean.unfed,
                                     mean.partial+mean.fully+mean.unfed,
                                     mean.partial+mean.fully+mean.unfed,
                                     mean.partial+mean.unfed),
          col=alpha("green", .5))
  
  
}


for(i in 1:length(treated.clusters)){
  unfed = bloodmeal_intervention %>%
    filter(Cluster == treated.clusters[i]) %>%
    mutate(propn_unfed = total_empty / total_evaluated) %>%
    pull(propn_unfed)
  
  partial = bloodmeal_intervention %>%
    filter(Cluster == treated.clusters[i]) %>%
    mutate(propn_partfed = total_half / total_evaluated) %>%
    pull(propn_partfed)
  
  fully = bloodmeal_intervention %>%
    filter(Cluster == treated.clusters[i]) %>%
    mutate(propn_fullfed = total_full / total_evaluated) %>%
    pull(propn_fullfed)
  
  total_evaluated = bloodmeal_intervention %>%
    filter(Cluster == treated.clusters[i]) %>%
    pull(total_evaluated)
  
  un = total_evaluated*unfed
  part = total_evaluated*partial
  full = total_evaluated*fully
  
  polygon(x=c(i+13, i+13, i+13-.25, i+13-.25), c(0, unfed, unfed, 0), col="red")
  polygon(x=c(i+13, i+13, i+13-.25, i+13-.25), c(unfed, partial+unfed,partial+unfed, unfed),
          col="blue")
  polygon(x=c(i+13, i+13, i+13-.25, i+13-.25), c(partial+unfed, partial+fully+unfed,
                                                 partial+fully+unfed,
                                                 partial+unfed),
          col="green")
  
  # uncertainty
  alpha.vec = c(1+un, 1+part, 1+full)
  uncertainty = rdirichlet(1000, alpha.vec)
  
  lines(x=c(i+13-.125, i+13-.125), y=c(quantile(uncertainty[,1], 0.025),
                                       quantile(uncertainty[,1], 0.975)))
  lines(x=c(i+13-.125, i+13-.125), y=c(quantile(uncertainty[,2], 0.025)+unfed,
                                       quantile(uncertainty[,2], 0.975)+unfed))
  lines(x=c(i+13-.125, i+13-.125), y=c(quantile(uncertainty[,3], 0.025)+
                                         unfed+partial,
                                       quantile(uncertainty[,3], 0.975)+
                                         unfed+partial))
  
  # plot the mode-predicted proportions
  mean.unfed = mean(prob.unfed.t)
  mean.partial = mean(prob.partfed.t)
  mean.fully = mean(prob.fullyfed.t)
  polygon(x=c(i+13+.25, i+13+.25, i+13, i+13), y=c(0, mean.unfed, mean.unfed, 0), 
          col=alpha("red", .5))
  polygon(x=c(i+13+.25, i+13+.25, i+13, i+13), c(mean.unfed, mean.partial+
                                                   mean.unfed,mean.partial+mean.unfed,
                                                 mean.unfed),
          col=alpha("blue", .5))
  polygon(x=c(i+13+.25, i+13+.25, i+13, i+13), c(mean.partial+mean.unfed,
                                                 mean.partial+mean.fully+mean.unfed,
                                                 mean.partial+mean.fully+mean.unfed,
                                                 mean.partial+mean.unfed),
          col=alpha("green", .5))
  
  
}
abline(v=13.5, lwd=4)

dev.off()
#### Mosquito Blood-Feeding- Baricentric plot ####

png("~/Dropbox/DengueEntomologicalEffects/results/ManuscriptPlots/PPDUnfed.png",
    width=1000, height=1000)

par(mfrow=c(6,6), mar=c(0,0,0,0), oma=c(0,0,0,0), cex=2)
for(i in 1:length(untreated.clusters)){
  
  
  triplot(x = NULL, y = NULL, z = NULL,
          frame = TRUE, label =F,
          grid = FALSE, center = FALSE,  set.par = TRUE)

  unfed = bloodmeal_baseline %>%
    filter(Cluster == untreated.clusters[i]) %>%
    mutate(propn_unfed = total_empty / total_evaluated) %>%
    pull(propn_unfed)
  
  part_fed = bloodmeal_baseline %>%
    filter(Cluster == untreated.clusters[i]) %>%
    mutate(propn_partfed = total_half / total_evaluated) %>%
    pull(propn_partfed)
  
  fully_fed = bloodmeal_baseline %>%
    filter(Cluster == untreated.clusters[i]) %>%
    mutate(propn_fullfed = total_full / total_evaluated) %>%
    pull(propn_fullfed)
  
  tripoints(unfed, part_fed, fully_fed, pch=19, col="red")

  unfed.quantiles.c = quantile(prob.unfed.c, c(0.025, .5, .975))
  partfed.quantiles.c = quantile(prob.partfed.c, c(0.025, .5, .975))
  fullfed.quantiles.c = quantile(prob.fullfed.c, c(0.025, .5, .975))
  
  tripoints(unfed.quantiles.c[2], partfed.quantiles.c[2], fullfed.quantiles.c[2],
            pch=19)
  tripoints(prob.unfed.c, prob.partfed.c, prob.fullfed.c)
  
}

plot.new()
plot.new()
plot.new()
plot.new()
plot.new()

for(i in 1:length(treated.clusters)){
  
  if(i == 2){
    plot.new()
    plot.new()
    plot.new()
    plot.new()
    plot.new()
  }
  
  triplot(x = NULL, y = NULL, z = NULL,
          frame = TRUE, label = FALSE,
          grid = FALSE, center = FALSE,  set.par = FALSE)
  
  unfed = bloodmeal_baseline %>%
    filter(Cluster == treated.clusters[i]) %>%
    mutate(propn_unfed = total_empty / total_evaluated) %>%
    pull(propn_unfed)
  
  part_fed = bloodmeal_baseline %>%
    filter(Cluster == treated.clusters[i]) %>%
    mutate(propn_partfed = total_half / total_evaluated) %>%
    pull(propn_partfed)
  
  fully_fed = bloodmeal_baseline %>%
    filter(Cluster == treated.clusters[i]) %>%
    mutate(propn_fullfed = total_full / total_evaluated) %>%
    pull(propn_fullfed)
  
  tripoints(unfed, part_fed, fully_fed, pch=19, col="red")
  
  unfed.quantiles.c = quantile(prob.unfed.c, c(0.025, .5, .975))
  partfed.quantiles.c = quantile(prob.partfed.c, c(0.025, .5, .975))
  fullfed.quantiles.c = quantile(prob.fullfed.c, c(0.025, .5, .975))
  
  tripoints(unfed.quantiles.c[2], partfed.quantiles.c[2], fullfed.quantiles.c[2],
            pch=19)
  tripoints(prob.unfed.c, prob.partfed.c, prob.fullfed.c)
  
}
dev.off()

#### Blood-Feeding Intervention ####
png("~/Dropbox/DengueEntomologicalEffects/results/ManuscriptPlots/PPDUnfedIntervention.png",
    width=1000, height=1000)
 
par(mfrow=c(6,6), mar=c(0,0,0,0), oma=c(0,0,0,0), cex=2)
for(i in 1:length(untreated.clusters)){
  
  
  triplot(x = NULL, y = NULL, z = NULL,
          frame = TRUE, label = F,
          grid = FALSE, center = FALSE,  set.par = TRUE)
  
  unfed = bloodmeal_intervention %>%
    filter(Cluster == untreated.clusters[i]) %>%
    mutate(propn_unfed = total_empty / total_evaluated) %>%
    pull(propn_unfed)
  
  part_fed = bloodmeal_intervention %>%
    filter(Cluster == untreated.clusters[i]) %>%
    mutate(propn_partfed = total_half / total_evaluated) %>%
    pull(propn_partfed)
  
  fully_fed = bloodmeal_intervention %>%
    filter(Cluster == untreated.clusters[i]) %>%
    mutate(propn_fullfed = total_full / total_evaluated) %>%
    pull(propn_fullfed)
  
  tripoints(unfed, part_fed, fully_fed, pch=19, col="red")
  
  unfed.quantiles.c = quantile(prob.unfed.c, c(0.025, .5, .975))
  partfed.quantiles.c = quantile(prob.partfed.c, c(0.025, .5, .975))
  fullfed.quantiles.c = quantile(prob.fullfed.c, c(0.025, .5, .975))
  
  tripoints(unfed.quantiles.c[2], partfed.quantiles.c[2], fullfed.quantiles.c[2],
            pch=19)
  tripoints(prob.unfed.c, prob.partfed.c, prob.fullfed.c)
  
}

plot.new()
plot.new()
plot.new()
plot.new()
plot.new()


for(i in 1:length(treated.clusters)){
  
  if(i == 2){
    plot.new()
    plot.new()
    plot.new()
    plot.new()
    plot.new()
  }
  
  triplot(x = NULL, y = NULL, z = NULL, 
          frame = TRUE, label = FALSE,
          grid = FALSE, center = FALSE,  set.par = TRUE)
  
  unfed = bloodmeal_intervention %>%
    filter(Cluster == treated.clusters[i]) %>%
    mutate(propn_unfed = total_empty / total_evaluated) %>%
    pull(propn_unfed)
  
  part_fed = bloodmeal_intervention %>%
    filter(Cluster == treated.clusters[i]) %>%
    mutate(propn_partfed = total_half / total_evaluated) %>%
    pull(propn_partfed)
  
  fully_fed = bloodmeal_intervention %>%
    filter(Cluster == treated.clusters[i]) %>%
    mutate(propn_fullfed = total_full / total_evaluated) %>%
    pull(propn_fullfed)
  
  tripoints(unfed, part_fed, fully_fed, pch=19, col="red")
  
  unfed.quantiles.t = quantile(prob.unfed.t, c(0.025, .5, .975))
  partfed.quantiles.t= quantile(prob.partfed.t, c(0.025, .5, .975))
  fullfed.quantiles.t = quantile(prob.fullfed.t, c(0.025, .5, .975))
  
  tripoints(unfed.quantiles.t[2], partfed.quantiles.t[2], fullfed.quantiles.t[2],
            pch=19)
  tripoints(prob.unfed.t, prob.partfed.t, prob.fullfed.t)
  
}



dev.off()

# legend
par(mfrow=c(1,1))
plot(0, type="n", xlim=c(0,1))
legend("center", c("Trial Data", "Model Prediction"), pch=19, col=c("red", "black"))

#### PPD bloodfeeding modified ####

par(mfrow=c(1,1))
plot(0, type = "n", xlim=c(1,26), ylim=c(-.1,1), 
     main="Mosquito Bloodfeeding, Baseline", xlab="", 
     xaxs="i", yaxs="i", ylab="Proportion")
rect.colors = c(rep(c("lightgrey", "white"), 14))
for(i in 1:28){
  rect(i, -.5, i+1, 600, col=rect.colors[i],
       border=NA)
}
abline(v=14)

for(i in 1:length(untreated.clusters)){
  
  unfed = bloodmeal_baseline %>%
    filter(Cluster == untreated.clusters[i]) %>%
    mutate(propn_unfed = total_empty / total_evaluated) %>%
    pull(propn_unfed)
  
  part_fed = bloodmeal_baseline %>%
    filter(Cluster == untreated.clusters[i]) %>%
    mutate(propn_partfed = total_half / total_evaluated) %>%
    pull(propn_partfed)
  
  fully_fed = bloodmeal_baseline %>%
    filter(Cluster == untreated.clusters[i]) %>%
    mutate(propn_fullfed = total_full / total_evaluated) %>%
    pull(propn_fullfed)
  
  rect(i, 0, i+.5, unfed, col="darkred")
  rect(i, unfed, i+.5, part_fed+unfed, col="red")
  rect(i, unfed+part_fed, i+.5,
       part_fed+unfed+fully_fed, col="darksalmon")
  
  rect(i+.5, 0, i+1, prob_inside_unfed_control,
       col="darkgreen")
  rect(i+.5, prob_inside_unfed_control,
       i+1, prob_inside_partfed_control+prob_inside_unfed_control, col="green")
  rect(i+.5,
       prob_inside_unfed_control+prob_inside_partfed_control,
       i+1,
       prob_inside_partfed_control+prob_inside_unfed_control+prob_inside_fullyfed_control,
       col="lightgreen")
  
}

for(i in 1:length(treated.clusters)){
  j=i+13
  unfed = bloodmeal_baseline %>%
    filter(Cluster == treated.clusters[i]) %>%
    mutate(propn_unfed = total_empty / total_evaluated) %>%
    pull(propn_unfed)
  
  part_fed = bloodmeal_baseline %>%
    filter(Cluster == treated.clusters[i]) %>%
    mutate(propn_partfed = total_half / total_evaluated) %>%
    pull(propn_partfed)
  
  fully_fed = bloodmeal_baseline %>%
    filter(Cluster == treated.clusters[i]) %>%
    mutate(propn_fullfed = total_full / total_evaluated) %>%
    pull(propn_fullfed)
  
  rect(j, 0, j+.5, unfed, col="darkred")
  rect(j, unfed, j+.5, part_fed+unfed, col="red")
  rect(j, unfed+part_fed, j+.5,
       part_fed+unfed+fully_fed, col="darksalmon")
  
  rect(j+.5, 0, j+1, prob_inside_unfed_control,
       col="darkgreen")
  rect(j+.5, prob_inside_unfed_control,
       j+1, prob_inside_partfed_control+prob_inside_unfed_control, col="green")
  rect(j+.5,
       prob_inside_unfed_control+prob_inside_partfed_control,
       j+1,
       prob_inside_partfed_control+prob_inside_unfed_control+prob_inside_fullyfed_control,
       col="lightgreen")
  
}

#### PPD bloodfeeding intervention ####
png("~/Dropbox/DengueEntomologicalEffects/results/ManuscriptPlots/PPDBloodfeeding.png")



par(mfrow=c(1,1))
plot(0, type = "n", xlim=c(1,26), ylim=c(-.1,1), 
     main="Mosquito Bloodfeeding, Intervention", xlab="", 
     xaxs="i", yaxs="i", ylab="Proportion")
rect.colors = c(rep(c("lightgrey", "white"), 14))
for(i in 1:28){
  rect(i, -.5, i+1, 600, col=rect.colors[i],
       border=NA)
}
abline(v=14)

for(i in 1:length(untreated.clusters)){
  
  unfed = bloodmeal_intervention %>%
    filter(Cluster == untreated.clusters[i]) %>%
    mutate(propn_unfed = total_empty / total_evaluated) %>%
    pull(propn_unfed)
  
  part_fed = bloodmeal_intervention %>%
    filter(Cluster == untreated.clusters[i]) %>%
    mutate(propn_partfed = total_half / total_evaluated) %>%
    pull(propn_partfed)
  
  fully_fed = bloodmeal_intervention %>%
    filter(Cluster == untreated.clusters[i]) %>%
    mutate(propn_fullfed = total_full / total_evaluated) %>%
    pull(propn_fullfed)
  
  rect(i, 0, i+.5, unfed, col="darkred")
  rect(i, unfed, i+.5, part_fed+unfed, col="red")
  rect(i, unfed+part_fed, i+.5,
       part_fed+unfed+fully_fed, col="darksalmon")
  
  rect(i+.5, 0, i+1, prob_inside_unfed_control,
       col="darkgreen")
  rect(i+.5, prob_inside_unfed_control,
       i+1, prob_inside_partfed_control+prob_inside_unfed_control, col="green")
  rect(i+.5,
       prob_inside_unfed_control+prob_inside_partfed_control,
       i+1,
       prob_inside_partfed_control+prob_inside_unfed_control+prob_inside_fullyfed_control,
       col="lightgreen")
  
}

for(i in 1:length(treated.clusters)){
  j=i+13
  unfed = bloodmeal_intervention %>%
    filter(Cluster == treated.clusters[i]) %>%
    mutate(propn_unfed = total_empty / total_evaluated) %>%
    pull(propn_unfed)
  
  part_fed = bloodmeal_intervention %>%
    filter(Cluster == treated.clusters[i]) %>%
    mutate(propn_partfed = total_half / total_evaluated) %>%
    pull(propn_partfed)
  
  fully_fed = bloodmeal_intervention %>%
    filter(Cluster == treated.clusters[i]) %>%
    mutate(propn_fullfed = total_full / total_evaluated) %>%
    pull(propn_fullfed)
  
  rect(j, 0, j+.5, unfed, col="darkred")
  rect(j, unfed, j+.5, part_fed+unfed, col="red")
  rect(j, unfed+part_fed, j+.5,
       part_fed+unfed+fully_fed, col="darksalmon")
  
  rect(j+.5, 0, j+1, prob_inside_unfed_treated,
       col="darkgreen")
  rect(j+.5, prob_inside_unfed_treated,
       j+1, prob_inside_partfed_treated+prob_inside_unfed_treated, col="green")
  rect(j+.5,
       prob_inside_unfed_treated+prob_inside_partfed_treated,
       j+1,
       prob_inside_partfed_treated+prob_inside_unfed_treated+
         prob_inside_fullyfed_treated,
       col="lightgreen")
  
}
dev.off()


prev.treated = samples[,paste0("X_", treated.clusters)]
prev.untreated = samples[,paste0("X_", (1:26)[-treated.clusters])]
median(colMeans(prev.treated)) / median(colMeans(prev.untreated))
par(mfrow=c(2,1))
hist(colMeans(prev.treated), xlim=c(0,.4))
hist(colMeans(prev.untreated), xlim=c(0,.4))

apply(prev.treated, 2, median)
apply(prev.treated, 2, median)