library(dplyr)

#### Does it make sense for epidemiology?

setwd("~/Dropbox/DengueEntomologicalEffects/code")
source("functions.R")

# cluster 1
epi=read.csv("epi.csv")
head(epi)
epi_1 = epi %>%
  filter(Cluster == 1)
foi_1 = epi %>%
  summarize(foi = sum(outcome) / sum(followup_time))
C=.45
# load prior particles
setwd("~/Dropbox/DengueEntomologicalEffects/results/Priors")
prior.samps = head(prior.samps, 1000)
foi_1_prior = rep(NA, 1000)
pi.U = rep(NA, 1000)
pi.T = rep(NA, 1000)
pi.tau = rep(NA, 1000)
for(i in 1:1000){
  mosquito_location_probs_treated = with(prior.samps[i,], location_probability(af, ap,
                                                                               qu, tau, d, rho, 
                                                         qt, alpha, C))
  pi.U = sum(mosquito_location_probs_treated[1:3])
  pi.tau = sum(mosquito_location_probs_treated[4:6])
  pi.T = sum(mosquito_location_probs_treated[7:9])
  gc_calc_i = with(prior.samps[i,],gc_calc(mu, qu, qt, gu, gt, gtau, tau, rho, C,
                                            pi.tau, pi.U, pi.T))
  
  ac_calc_i = with(prior.samps[i,], ac_calc(af, ap, rho, C, alpha, qt, qu, tau))
  
  foi_1_prior[i] = with(samples[i,], force_of_infection(b, lambda, ac_calc_i, c, X_1, gc_calc_i, n))
}

# find the likelihood of each particles

#### Likelihood Function ####

likelihood_epi_1=function(params){
  
  # SR-specific parameters
  alpha = params[1]
  qt = params[2]
  gt=params[3]
  mu=params[4]
  rho=params[5]
  
  # mosquito parameters
  lambda = params[6]
  af = params[7]
  ap = params[8]
  qu = params[9]
  c = params[10]
  gu=params[11]
  n=params[12]
  tau=params[13]
  
  # human parameters
  b=params[14]
  X=params[15:40]
  
  gtau = params[41]
  phi=params[42]
  
  d=params[43]
  
  catch_prop = params[44]
  
  # biting rate, and modified biting rates
  a = af+ap
  ac = ac_calc(af, ap, rho, C, alpha, qt, qu, tau)
  ac_control = ac_calc(af, ap, rho, C=0, alpha, qt, qu, tau)
  
  # find the probability of mosquitoes being located in each compartment, control
  # and treatment
  mosquito_location_probs_control = location_probability(af, ap, qu, tau, d, rho, 
                                                         qt, alpha, C=0)[1:6]
  mosquito_location_probs_treated = location_probability(af, ap, qu, tau, d, rho, 
                                                         qt, alpha, C)
  
  # overall inside a treated house probability and 
  # conditional blood meal probabilities, subbing these into the vectors
  # only in treated clusters
  prob_inside_treated =  mosquito_location_probs_treated[7]+ 
    mosquito_location_probs_treated[8]+ 
    mosquito_location_probs_treated[9] 
  
  prob_transition_treated = mosquito_location_probs_treated[4]+
    mosquito_location_probs_treated[5]+
    mosquito_location_probs_treated[6]
  
  # probability of being in an untreated house in treatment cluster
  prob_inside_untreated_house_treated_cluster =  mosquito_location_probs_treated[1]+ 
    mosquito_location_probs_treated[2]+ 
    mosquito_location_probs_treated[3] 
  
  
  ## probability of death, and modified death
  gc = gc_calc(mu, qu, qt, gu, gt,tau, rho, C, gtau, 
               pi.tau=prob_transition_control, pi.U=prob_inside_control, pi.T=0)
  gc_control = gc_calc(mu=0, qu, qt=qu, gu, gt=gu,tau, rho, C, gtau,
                       pi.tau=prob_transition_treated,
                       pi.U = prob_inside_untreated_house_treated_cluster,
                       pi.T=prob_inside_treated)
  
  # force of infection under 0 coverage and under C
  foi = force_of_infection(b, lambda, ac_control, c, X, gc_control, n)
  foi_c = force_of_infection(b, lambda, ac, c, X, gc, n) 
  
  foi_trial = foi
  foi_trial[treated.clusters] = foi_c[treated.clusters]
  
  # parity rate
  parity_rate = ac_control / (ac_control+gc_control)
  parity_rate_c = rep(parity_rate, 26)
  parity_rate_c[treated.clusters] = ac / (ac + gc)
  
  
  # mosquito to human ratio in untreated and treated clusters
  m_control=lambda/gc_control
  mc_trt = lambda/gc
  
  # number of mosquitoes per household times probability of being inside 
  # (this is the expected number of mosquitoes one would expect to find in a 
  # 100% efficient aspiration indoors), control and treatment
  f_control = p * m_control * prob_inside_control
  fc_trt = p * mc_trt * prob_inside_treated
  
  # likelihood for incidence data
  # flip the 0's and 1's since we are looking at the probability of not seroconverting
  epi_likelihood = sum(dbinom(1-epi$outcome, size=rep(1, nrow(epi)), 
                              prob= exp(-foi_trial[epi$Cluster]*epi$followup_time),
                              log=T))
  
  # likelihood for baseline parity
  
  baseline_parity_likelihood = sum(dbinom(parity_baseline$parous, 
                                          parity_baseline$total_caught,
                                          prob=parity_rate, log=T))
  intervention_parity_likelihood = sum(dbinom(parity_intervention$parous, 
                                              parity_intervention$total_caught,
                                              prob=parity_rate_c, log=T))
  
  # likelihood for mosquito density
  baseline_abundance_likelihood = sum(dnbinom(abundance_baseline$total_caught,
                                              mu=f_control*catch_prop, size=phi, log=T))
  intervention_abundance_control_likelihood = sum(
    dnbinom(abundance_intervention_control$total_caught,
            mu=f_control*catch_prop,
            size=phi, log=T))
  intervention_abundance_treatment_likelihood = sum(
    dnbinom(abundance_intervention_treatment$total_caught,
            mu=fc_trt*catch_prop, size=phi, log=T))
  
  # likelihood for blood fed status
  
  baseline_fullyfed_likelihood = sum(dbinom(bloodmeal_baseline$total_full,
                                            size=bloodmeal_baseline$total_evaluated,
                                            prob=prob_inside_fullyfed_control,
                                            log=TRUE))
  
  
  
  baseline_partfed_likelihood = sum(dbinom(bloodmeal_baseline$total_half,
                                           size=bloodmeal_baseline$total_evaluated,
                                           prob=prob_inside_partfed_control,
                                           log=TRUE))
  baseline_unfed_likelihood = sum(dbinom(bloodmeal_baseline$total_empty,
                                         size=bloodmeal_baseline$total_evaluated,
                                         prob=prob_inside_unfed_control,
                                         log=TRUE))
  
  intervention_fullyfed_likelihood = sum(dbinom(bloodmeal_intervention$total_full,
                                                size=bloodmeal_intervention$total_evaluated,
                                                prob=prob_inside_fullyfed_treated[bloodmeal_intervention$Cluster],
                                                log=TRUE))
  intervention_partfed_likelihood = sum(dbinom(bloodmeal_intervention$total_half,
                                               size=bloodmeal_intervention$total_evaluated,
                                               prob=prob_inside_partfed_treated[bloodmeal_intervention$Cluster],
                                               log=TRUE))
  intervention_unfed_likelihood = sum(dbinom(bloodmeal_intervention$total_empty,
                                             size=bloodmeal_intervention$total_evaluated,
                                             prob=prob_inside_unfed_treated[bloodmeal_intervention$Cluster],
                                             log=TRUE))
  
  return(epi_likelihood + baseline_parity_likelihood + intervention_parity_likelihood +
           intervention_abundance_control_likelihood + baseline_abundance_likelihood+
           intervention_abundance_treatment_likelihood + 
           baseline_fullyfed_likelihood + baseline_partfed_likelihood +
           baseline_unfed_likelihood + intervention_fullyfed_likelihood + 
           intervention_partfed_likelihood +
           intervention_unfed_likelihood)
  
}
