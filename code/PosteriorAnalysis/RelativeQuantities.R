#### import samples and functions ####

setwd("~/Dropbox/DengueEntomologicalEffects/results/Posteriors")
load("~/Dropbox/DengueEntomologicalEffects/results/Posteriors/SequentialMonteCarlo.RData")

setwd("~/Dropbox/DengueEntomologicalEffects/code")
source("functions.R")

# priors 
setwd("~/Dropbox/DengueEntomologicalEffects/results/Priors")
load("~/Dropbox/DengueEntomologicalEffects/results/Priors/PriorSamples.RData")

#### compute ac and gc ####

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

C= 0.45

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


modified_biting = with(samples, ac_calc(af, ap, rho, 0.45, alpha, qt, qu, tau))
baseline_biting = with(samples, ac_calc(af, ap, rho=0, 0, alpha=1, qu, qu, tau))
relative_biting = modified_biting / baseline_biting

modified_death = with(samples, gc_calc(mu, qu, qt, gu, gt, gtau, tau, rho, .45, pi.tau.treated, pi.U.treated, pi.T.treated))
baseline_death = with(samples, gc_calc(mu, qu, qu, gu, gu, gtau, tau, rho=0, C, pi.tau.control, pi.U.treated, 0))
relative_death = modified_death / baseline_death

# compute on priors
modified_biting_prior = with(prior.samps, ac_calc(af, ap, rho, 0.45, alpha, qt,
                                                  qu, tau))
baseline_biting_prior = with(prior.samps, ac_calc(af, ap, rho=0, 0, alpha=1, qu,
                                                  qu, tau))
relative_biting_prior = modified_biting_prior / baseline_biting_prior

modified_death_prior = with(prior.samps, gc_calc(mu, qu, qt, gu, gt, gtau, tau, rho, .45))
baseline_death_prior = with(prior.samps, gc_calc(mu, qu, qu, gu, gu, gtau, tau, rho=0, C))
relative_death_prior = modified_death_prior / baseline_death_prior

png("~/Dropbox/DengueEntomologicalEffects/results/ManuscriptPlots/OverallEffects.png",
    width=1000, height=600)

par(mfrow=c(1,2), cex=1.2)
plot(density(relative_biting), col="blue", main="Relative Biting Rate",
     xlim=c(0.4,1.5), xlab="")
polygon(density(relative_biting), col="blue")
polygon(density(relative_biting_prior), col="red")
abline(v=1)
legend("topright", c("Prior", "Posterior"), col=c("red", "blue"), pch=19)

plot(density(relative_death), col="blue", main="Relative Death Rate",
     xlim=c(0.8,1.2), xlab="")
polygon(density(relative_death), col="blue")
polygon(density(relative_death_prior), col="red")
abline(v=1)
dev.off()

quantile(1-relative_biting, c(0.025, .5, .975))
quantile(1-relative_death, c(0.025, .5, .975))

mean(with(samples, gc_calc(mu, qu, qt, gu, .4, gtau, tau, rho, .45)))
mean(with(samples, gc_calc(mu, qu, qt, gu, .5, gtau, tau, rho, .45)))

#### Relative FOI, parity, expected catch number ####

baseline_foi= with(samples, force_of_infection(b, lambda, af+ap, c, X_1, gu, n))
intervention_foi = with(samples, force_of_infection(b, lambda, modified_biting, c, X_1, modified_death, n))
relative_foi = intervention_foi / baseline_foi

parity_foi = with(samples, parity(af+ap, gu))

par(mfrow=c(2,2))
plot(density(relative_foi), col="blue", "Relative FOI", xlab="")
polygon(density(relative_foi), col="blue")



