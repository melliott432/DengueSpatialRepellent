library(tidyverse)

#### import modules ####

load("~/Dropbox/DengueEntomologicalEffects/results/Posteriors/SequentialMonteCarlo.RData")

load("~/Dropbox/DengueEntomologicalEffects/results/Priors/PriorSamples.RData")

prior.samps = prior.samps %>%
  dplyr::select(c(alpha,qt,gt, mu, rho, lambda, af,ap,
                  qu , c, gu,n,tau, b,paste0("X_", 1:26)), gtau, phi, d, catch_prop)

samples=head(samples,1000)
samples$time_untreated = 1 / samples$qu
samples$time_treated = 1/ samples$qt
samples$relative_exit_rate = samples$qt / samples$qu
samples$relative_mortality_rate = samples$gt / samples$gu
samples$time_spent_in_transit = samples$tau
samples$baseline_lifespan = 1 / samples$gu

prior.samps$time_untreated = 1 / prior.samps$qu
prior.samps$time_treated = 1/ prior.samps$qt
prior.samps$relative_exit_rate = prior.samps$qt / prior.samps$qu

prior.samps$relative_mortality_rate = prior.samps$gt / prior.samps$gu
prior.samps$time_spent_in_transit = prior.samps$tau
prior.samps$baseline_lifespan = 1 / prior.samps$gu
prior.samps$mortality_during_transition = 1 / prior.samps$gtau

#### Divide columns, make names ####

mains = c("Relative Feeding Rate", "Exit Rate- Treated", "Mortality-treated",
          "Probability of immediate death",
          "Repellency Probability", "Mosquito Emergence", "Full Blood-Feeding Rate",
          "Partial Blood-Feeding Rate", "Exit Rate- Untreated",
          "Human-to-Mosquito Transmission Probability",
          "Baseline Mortality", "Entomological Incubation Period",
          "Transition probability", "Mosquito-to-Human Transmission Probability",
          "Human Recovery Rate", paste0("Prevalence in Cluster ", 1:26),
          "Mortality Rate during Transition",
          "Dispersion Parameter for Mosquito Abundance",
          "Time spent in untreated clusters",
          "Time spent in treated clusters", "Relative exit rate in treated clusters",
          "Relative mortality rate in treated clusters", "Infectious period", 
          "Time spent in transit", "Baseline Lifespan", "Mortality During Transition")

# divide parameters into groups for density plots

#### SR related parameters ####


# xlimit values
xlim.bottom=rep
xlim.top=c(2, .7, 1.5, 1, 1)
ylim.top=c(4, 8, 8, 4, 5)

#### SR-related parameters ####

SR_related_params = samples %>%
  dplyr::select(alpha, mu, rho, relative_exit_rate, relative_mortality_rate)
SR_related_priors = prior.samps %>%
  dplyr::select(alpha, mu, rho, relative_exit_rate, relative_mortality_rate)
mains=c("Relative Feeding Rate in \n treated homes", "Probability of \n immediate death",
        "Repellency Probability", "Relative exit rate\n in treated homes", 
        "Relative mortality rate\n in treated homes")
xlim.bottom = c(0,0.05,0,0.1,0)
xlim.top=c(2, .2, 1, 5, 5)
ylim.top = c(5, 50, 15, 4 ,2)
vline = c(1, 0, 0, 1, 1)

png(file="~/Dropbox/DengueEntomologicalEffects/results/ManuscriptPlots/SRparams.png",
    height=700, width=800)
par(mfrow=c(2,3), mar=c(2, 2, 3, 1), cex=1.3)
for(i in 1:ncol(SR_related_params)){
  plot(density(SR_related_params[,i]), col="blue", main=mains[i], xlab="Values",
       xlim=c(xlim.bottom[i], xlim.top[i]), ylim=c(0, ylim.top[i]), cex=2)
  polygon(density(SR_related_priors[,i]), col="red")
  polygon(density(SR_related_params[,i]), col="blue")
  abline(v=vline[i])
}
legend("topright", c("Prior", "Posterior"), col=c("red", "blue"), pch=19)
dev.off()

# 95% CI
1-quantile(samples$alpha, c(.025, .5, .975))
1-quantile(samples$relative_exit_rate, c(.025, .5, .975))
quantile(samples$mu, c(.025, .5, .975))
quantile(samples$relative_mortality_rate, c(.025, .5, .975))
quantile(samples$rho,  c(.025, .5, .975))

#### Non-SR related parameters #####

samples$time_spent_untreated = 1 / samples$qu
samples$digestion_time = 1/samples$d

prior.samps$time_spent_untreated = 1 / prior.samps$qu
prior.samps$digestion_time = 1/prior.samps$d

png(file="~/Dropbox/DengueEntomologicalEffects/results/ManuscriptPlots/NonSRparams.png",
    height=1000, width=1000)
par(mfrow=c(3,3), cex=1.2, mar=c(3,3,3,0))

# mosquito emergence
plot(density(samples$lambda), col="blue", main="Mosquito \nEmergence",
     xlab="Values", cex=2,
     xlim=c(0,10))
polygon(density(prior.samps$lambda), col="red")
polygon(density(samples$lambda), col="blue")

plot(density(samples$af), col="blue", main="Full Blood-Feeding Rate",
     xlab="Values", cex=2,
     xlim=c(0.01,0.2))
polygon(density(prior.samps$af), col="red")
polygon(density(samples$af), col="blue")

plot(density(samples$ap), col="blue", main="Partial Blood-Feeding Rate",
     xlab="Values", cex=2,
     xlim=c(0,0.3))
polygon(density(prior.samps$ap), col="red")
polygon(density(samples$ap), col="blue")

plot(density(samples$time_untreated), col="blue",
     main="Time spent in untreated\n houses (days)",
     xlab="Values", cex=2,
     xlim=c(1,8))
polygon(density(prior.samps$time_untreated), col="red")
polygon(density(samples$time_untreated), col="blue")

plot(density(samples$n), col="blue",
     main="Entomological \nIncubation \n Period",
     xlab="Values", cex=2,
     xlim=c(4,11))
polygon(density(prior.samps$n), col="red")
polygon(density(samples$n), col="blue")

plot(density(samples$c), col="blue",
     main="Human to Mosquito \nTransmission Probability",
     xlab="Values", cex=2,
     xlim=c(0.3, 1))
polygon(density(prior.samps$c), col="red")
polygon(density(samples$c), col="blue")

plot(density(samples$b), col="blue",
     main="Mosquito to Human \n Transmission Probability",
     xlab="Values", cex=2,
     xlim=c(0, 1))
polygon(density(prior.samps$b), col="red")
polygon(density(samples$b), col="blue")

plot(density(samples$digestion_time), col="blue",
     main="Digestion Time for Mosquitoes \n (Days)",
     xlab="Values", cex=2,
     xlim=c(3,6))
polygon(density(prior.samps$digestion_time), col="red")
polygon(density(samples$digestion_time), col="blue")

plot(density(samples$catch_prop), col="blue",
     main="Aspiration Catch Efficiency",
     xlab="Values", cex=2,
     xlim=c(0.01, .08))
polygon(density(prior.samps$catch_prop), col="red")
polygon(density(samples$catch_prop), col="blue")


dev.off()

##### Prevalences by Cluster ####


prevalence_params = samples %>%
  dplyr::select(contains("X")) %>%
  dplyr::select(-relative_exit_rate)
prevalence_priors = prior.samps %>%
  dplyr::select(contains("X"))  %>%
  dplyr::select(-relative_exit_rate)
mains=paste("Prevalence in \n Cluster ", 1:26)

# xlimit values

xlim.bottom=rep(0, 26)
xlim.top=c(0.6, 0.2, 0.5, 0.4)
ylim.top = rep(100, 26)
png(file="~/Dropbox/DengueEntomologicalEffects/results/ManuscriptPlots/Prevalenceparams.png",
    height=700, width=700)
par(mfrow=c(4,4), mar=c(2,2,2,2))
for(i in 1:16){
  plot(density(samples[,paste0("X_", i)]), col="blue", main="", xlab="Values",
       xlim=c(0, 1), ylim=c(0, 40), cex=2)
  polygon(density(prior.samps[,paste0("X_", i)]), col="red")
  polygon(density(samples[,paste0("X_", i)]), col="blue")
  
}
legend("topright", c("Prior", "Post"), col=c("red", "blue"), pch=19)
dev.off()

#### mosquito locations ####

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


#### modified blood-feeding and mortality rates ####

setwd("~/Dropbox/DengueEntomologicalEffects/code")
source("functions.R")

# compute
ac_calc_control = with(samples, ac_calc(af, ap, rho, C=0, alpha, qt, qu, tau))
ac_calc_treatment = with(samples, ac_calc(af, ap, rho, C=.42, alpha, qt, qu, tau))
relative_biting_rate = ac_calc_treatment/ ac_calc_control

gc_calc_control = with(samples, gc_calc(mu, qu, qt, gu, gt, gtau, tau, rho, C=0, pi.tau.control,
                       pi.U.control, pi.T.control))
gc_calc_treatment = with(samples, gc_calc(mu, qu, qt, gu, gt, gtau, tau, rho,
                                        C=.42, pi.tau.treated,
                                        pi.U.treated, pi.T.treated))
relative_mortality_rate = gc_calc_treatment/ gc_calc_control
png(file="~/Dropbox/DengueEntomologicalEffects/results/ManuscriptPlots/ModifiedValues.png",
    height=500, width=1000)
par(mfrow=c(1,2), cex=1.2)
plot(density(relative_biting_rate), xlim=c(0.85, 1.05),
     col="blue",main="Relative biting rate \nin treated clusters",
     xlab="")
polygon(density(relative_biting_rate), col="blue")
abline(v=1)
plot(density(relative_mortality_rate), col="blue",xlim=c(1.4, 2),
     main="Relative mosquito mortality \nin treated clusters",
     xlab="")
polygon(density(relative_mortality_rate), col="blue")
abline(v=1)
dev.off()

1-quantile(relative_biting_rate, c(0.025, .5, 0.975))
quantile(relative_mortality_rate, c(0.025, .5, 0.975))
