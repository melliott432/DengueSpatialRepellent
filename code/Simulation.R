
library(tidyverse)
library(bbmle)
library(doParallel)

#### load data ####

epi = read.csv(file="epi.csv")
parity_baseline = read.csv(file="parity_baseline.csv")
parity_intervention = read.csv(file="parity_intervention.csv")
load("treatedclusters.RData")
untreated.clusters=(1:26)[-treated.clusters]
abundance_baseline = read.csv(file="abundance_baseline.csv")
abundance_intervention_control = read.csv(file="abundance_intervention_control.csv")
abundance_intervention_treatment = read.csv(file="abundance_intervention_treatment.csv")
bloodmeal_baseline = read.csv(file="bloodmeal_baseline.csv")
bloodmeal_intervention = read.csv(file="bloodmeal_intervention.csv")

C = .45
p=5.2

source("functions.R")







#### Priors, upper and lower ####

load("PriorSamples.RData")
prior.samps = prior.samps %>%
  # mutate(phi = log(phi)) %>%
  dplyr::select(c(alpha,qt,gt,mu,rho,lambda,af, ap,qu ,c,gu,n,tau,b,paste0("X_", 1:26), gtau,
                phi, d, catch_prop))
upper=apply(prior.samps, 2, max)
upper["ap"] = 2
upper["af"] = 1
lower=apply(prior.samps, 2, min)
lower["lambda"]=.001
lower["catch_prop"] = .00001
lower["phi"] = .0001
lower["gu"] = lower["gt"]
upper["gu"] = upper["gt"]
lower["gt"] = lower["gu"]
upper["qt"] = upper["qu"]
lower[15:40] = .00000001
upper[15:40] = .2

prior.samps=prior.samps
for(i in 1:nrow(prior.samps)){
  prior.samps[i,] = runif(44, lower, upper)
}


#### Simulate the data sets ####


simulate_data = function(X_1, X_2, X_3, X_4, X_5, X_6, X_7, X_8, X_9, X_10,
                         X_11, X_12, X_13, X_14, X_15, X_16, X_17, X_18, X_19, X_20,
                         X_21, X_22, X_23, X_24, X_25, X_26,af,ap,gu,gtau,
                         qu, tau,n, b,c,alpha, mu,
                         gt,rho,qt,lambda,d,
                         catch_prop,phi){
  X=c(X_1, X_2, X_3, X_4, X_5, X_6, X_7, X_8, X_9, X_10,
      X_11, X_12, X_13, X_14, X_15, X_16, X_17, X_18, X_19, X_20,
      X_21, X_22, X_23, X_24, X_25, X_26)
  
  
  
  
  epi_sim = epi
  
  bloodmeal_baseline_sim = bloodmeal_baseline
  bloodmeal_intervention_sim = bloodmeal_intervention
  
  parity_baseline_sim = parity_baseline
  parity_intervention_sim = parity_intervention
  
  abundance_baseline_sim = abundance_baseline
  abundance_baseline_sim$total_caught = NA
  
  abundance_intervention_control_sim = abundance_intervention_control
  abundance_intervention_control_sim$total_caught = NA
  
  abundance_intervention_treatment_sim = abundance_intervention_treatment
  abundance_intervention_treatment_sim$total_caught = NA
  
  #### Quantities
  
  
  ac_control = ac_calc(af, ap, rho, C=0, alpha, qt, qu, tau)
  ac_treatment = ac_calc(af, ap, rho, C=.45, alpha, qt, qu, tau)
  
  # mosquito locations
  locations_control = location_probability(af, ap, qu, tau, d, rho, qt, alpha, C=0)
  locations_treatment = location_probability(af, ap, qu, tau, d, rho, qt, alpha, C)
  
  prob_inside_control = sum(locations_control[1:3])
  prob_trans_control = sum(locations_control[4:6])
  
  prob_inside_untreated_treatment = sum(locations_treatment[1:3])
  prob_trans_treatment = sum(locations_treatment[4:6])
  prob_inside_treated_treatment = sum(locations_treatment[7:9])
  
  # blood feeding status
  prob_fullyfed_control = locations_control[1]/prob_inside_control
  prob_partfed_control = locations_control[2]/prob_inside_control
  prob_unfed_control = locations_control[3]/prob_inside_control
  
  prob_fullyfed_treatment = locations_treatment[7]/prob_inside_treated_treatment
  prob_partfed_treatment = locations_treatment[8]/prob_inside_treated_treatment
  prob_unfed_treatment = locations_treatment[9]/prob_inside_treated_treatment
  
  # mortality
  gc_control = gc_calc(mu, qu, qt, gu, gt, gtau, tau, rho, C=0, pi.tau=prob_trans_control,
                       pi.U = prob_inside_control, pi.T=0)
  gc_treatment = gc_calc(mu, qu, qt, gu, gt, gtau, tau, rho, C, pi.tau=prob_trans_treatment,
                         pi.U = prob_inside_untreated_treatment, pi.T=prob_inside_treated_treatment)
  
  foi_control = force_of_infection(b, lambda, ac_control, c, X, gc_control, n)
  foi_treatment = force_of_infection(b, lambda, ac_treatment, c, X, gc_treatment, n)
  foi_trial = foi_control
  foi_trial[treated.clusters] = foi_treatment[treated.clusters]
  
  # parity 
  parity_control = parity(ac_control, gc_control)
  parity_treatment = parity(ac_treatment, gc_treatment)
  
  # ratio of mosquitoes to humans
  m_control = lambda/gc_control
  m_treatment = lambda / gc_treatment
  
  #### Generate epi data
  # probability of not seroconverting
  prob_epi= exp(-foi_trial[epi_sim$Cluster]*epi_sim$followup_time)
  # convert to probability of seroconverting
  prob_epi = 1-prob_epi
  epi_sim$outcome = rbinom(nrow(epi_sim), size=1, 
                           prob= prob_epi)
  
  #### Generate bloodmeal data ####
  for(i in 1:26){
    bloodmeal_baseline_sim[i,2:4] = rmultinom(1, size=bloodmeal_baseline_sim$total_evaluated[i],
                                              prob=c(prob_fullyfed_control, prob_partfed_control,
                                                     prob_unfed_control))
  }
  
  #### Intervention Period
  
  
  ## Control Arm
  for(i in untreated.clusters){
    bloodmeal_intervention_sim[i,2:4] = rmultinom(1, size=bloodmeal_intervention_sim$total_evaluated[i],
                                                  prob=c(prob_fullyfed_control, prob_partfed_control,
                                                         prob_unfed_control))
  }
  ## Treatment Arm
  
  for(i in treated.clusters){
    bloodmeal_intervention_sim[i,2:4] = rmultinom(1, size=bloodmeal_intervention_sim$total_evaluated[i],
                                                  prob=c(prob_fullyfed_treatment, prob_partfed_treatment,
                                                         prob_unfed_treatment))
  }
  
  #### parity data ####
  
  parity_baseline_sim$nulliparous =parity_baseline_sim$parous = NA
  parity_baseline_sim$parous = rbinom(nrow(parity_baseline), size=parity_baseline_sim$total_caught,
                                      prob=parity_control)
  
  parity_intervention_sim$nulliparous =parity_intervention_sim$parous = NA
  parity_intervention_sim$parous[untreated.clusters] = rbinom(nrow(parity_intervention[untreated.clusters,]),
                                                              size=parity_intervention_sim[untreated.clusters,]$total_caught,
                                                              prob=parity_control)
  
  ## Treatment
  parity_intervention_sim$parous[treated.clusters] = rbinom(nrow(parity_intervention[treated.clusters,]),
                                                            size=parity_intervention_sim[treated.clusters,]$total_caught,
                                                            prob=parity_treatment)
  # abundance
  
  abundance_baseline_sim$total_caught = rnbinom(nrow(abundance_baseline_sim),
                                                mu=prob_inside_control*p*m_control*catch_prop,
                                                size=phi)
  abundance_intervention_control_sim$total_caught = rnbinom(nrow(abundance_intervention_control_sim),
                                                            mu=prob_inside_control*p*m_control*catch_prop,
                                                            size=phi)
  abundance_intervention_treatment_sim$total_caught = rnbinom(nrow(abundance_intervention_treatment_sim),
                                                              mu=prob_inside_treated_treatment*p*
                                                                m_treatment*catch_prop,
                                                              size=phi)
  
  return(list(abundance_baseline_sim, abundance_intervention_control_sim, 
              abundance_intervention_treatment_sim, parity_baseline_sim, parity_intervention_sim,
              bloodmeal_baseline_sim, bloodmeal_intervention_sim, epi_sim))
}

#### use this data set to fidn the likelihood ####

# want to return: true and inferred values of the parameters
find_max_likelihood = function(X_1, X_2, X_3, X_4, X_5, X_6, X_7, X_8, X_9, X_10,
                               X_11, X_12, X_13, X_14, X_15, X_16, X_17, X_18, X_19, X_20,
                               X_21, X_22, X_23, X_24, X_25, X_26,af,ap,gu,gtau,
                               qu, tau,n, b,c,alpha, mu,
                               gt,rho,qt,lambda,d,
                               catch_prop,phi){
  
   simulated_data = simulate_data(X_1, X_2, X_3, X_4, X_5, X_6, X_7, X_8, X_9, X_10,
                                X_11, X_12, X_13, X_14, X_15, X_16, X_17, X_18, X_19, X_20,
                                X_21, X_22, X_23, X_24, X_25, X_26,af,ap,gu,gtau,
                                qu, tau,n, b,c,alpha, mu,
                                gt,rho,qt,lambda,d,
                                catch_prop,phi)
   
   abundance_baseline_sim = simulated_data[[1]]
   abundance_intervention_control_sim = simulated_data[[2]]
   abundance_intervention_treatment_sim = simulated_data[[3]]
   parity_baseline_sim= simulated_data[[4]]
   parity_intervention_sim= simulated_data[[5]]
   bloodmeal_baseline_sim = simulated_data[[6]]
   bloodmeal_intervention_sim = simulated_data[[7]]
   epi_sim= simulated_data[[8]]
                               
                               
  
   
   likelihood_sim=function(alpha, qt, gt, mu, rho, lambda, af, ap, qu, c, gu, n, tau, b, 
                           X_1, X_2, X_3, X_4, X_5, X_6, X_7, X_8, X_9, X_10, X_11, X_12, 
                           X_13, X_14, X_15, X_16, X_17, X_18, X_19, X_20, X_21, 
                           X_22, X_23, X_24, X_25, X_26, gtau, phi, d, catch_prop){
     
     # SR-specific parameters
     X=c(X_1, X_2, X_3, X_4, X_5, X_6, X_7, X_8, X_9, X_10, X_11, X_12, 
         X_13, X_14, X_15, X_16, X_17, X_18, X_19, X_20, X_21, 
         X_22, X_23, X_24, X_25, X_26)
     
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
     
     # add up probabilities of being inside and find conditional bloodmeal
     # status probabilities
     prob_inside_control =  mosquito_location_probs_control[1]+ 
       mosquito_location_probs_control[2]+ 
       mosquito_location_probs_control[3]
     prob_inside_fullyfed_control = mosquito_location_probs_control[1] / prob_inside_control
     prob_inside_partfed_control = mosquito_location_probs_control[2] / prob_inside_control
     prob_inside_unfed_control = mosquito_location_probs_control[3] / prob_inside_control
     
     # duplicate the vectors before adding treatment condition values for 
     # treated clusters
     # note control values are the same for every cluster
     prob_inside_fullyfed_treated = prob_inside_fullyfed_control = 
       rep(prob_inside_fullyfed_control, 26)
     prob_inside_partfed_treated = prob_inside_partfed_control = 
       rep(prob_inside_partfed_control, 26)
     prob_inside_unfed_treated = prob_inside_unfed_control = 
       rep(prob_inside_unfed_control, 26)
     
     # overall inside a treated house probability and 
     # conditional blood meal probabilities, subbing these into the vectors
     # only in treated clusters
     prob_inside_treated =  mosquito_location_probs_treated[7]+ 
       mosquito_location_probs_treated[8]+ 
       mosquito_location_probs_treated[9] 
     if(prob_inside_treated == 0){
       prob_inside_fullyfed_treated[treated.clusters] = 0
       prob_inside_partfed_treated[treated.clusters] = 0
       prob_inside_unfed_treated[treated.clusters] =0
     } else{
       prob_inside_fullyfed_treated[treated.clusters] = mosquito_location_probs_treated[7]/ 
         prob_inside_treated
       prob_inside_partfed_treated[treated.clusters] = mosquito_location_probs_treated[8]/
         prob_inside_treated
       prob_inside_unfed_treated[treated.clusters] = mosquito_location_probs_treated[9]/
         prob_inside_treated
     }
     
     
     # calculate the probability that mosquitoes are in transition
     prob_transition_control = mosquito_location_probs_control[4]+
       mosquito_location_probs_control[5]+
       mosquito_location_probs_control[6]
     
     prob_transition_treated = mosquito_location_probs_treated[4]+
       mosquito_location_probs_treated[5]+
       mosquito_location_probs_treated[6]
     
     # probability of being in an untreated house in treatment cluster
     prob_inside_untreated_house_treated_cluster =  mosquito_location_probs_treated[1]+ 
       mosquito_location_probs_treated[2]+ 
       mosquito_location_probs_treated[3] 
     
     
     ## probability of death, and modified death
     gc = gc_calc(mu, qu, qt, gu, gt,tau, rho, C, gtau, 
                  pi.tau=prob_transition_treated,
                  pi.U=prob_inside_untreated_house_treated_cluster, pi.T=prob_inside_treated)
     gc_control = gc_calc(mu=0, qu, qt=qu, gu, gt=gu,tau, rho, C=0, gtau,
                          pi.tau=prob_transition_control,
                          pi.U = prob_inside_control,
                          pi.T=0)
     
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
     fc_trt[fc_trt==0] = 0.000001
     f_control[fc_trt==0] = 0.000001
     
     # likelihood for incidence data
     # flip the 0's and 1's since we are looking at the probability of not seroconverting
     prob_epi= exp(-foi_trial[epi$Cluster]*epi$followup_time)
     
     prob_epi[prob_epi == 1] = 1 - .0000001
     prob_epi[prob_epi == 0] = 0+ .0000001
     
     epi_likelihood = sum(dbinom(1-epi_sim$outcome, size=rep(1, nrow(epi)), 
                                 prob= prob_epi,
                                 log=T))
     
     # likelihood for baseline parity
     
     baseline_parity_likelihood = sum(dbinom(parity_baseline_sim$parous, 
                                             parity_baseline_sim$total_caught,
                                             prob=parity_rate, log=T))
     intervention_parity_likelihood = sum(dbinom(parity_intervention_sim$parous, 
                                                 parity_intervention_sim$total_caught,
                                                 prob=parity_rate_c, log=T))
     
     # likelihood for mosquito density
     baseline_abundance_likelihood = sum(dnbinom(abundance_baseline_sim$total_caught,
                                                 mu=f_control*catch_prop, size=phi, log=T))
     intervention_abundance_control_likelihood = sum(
       dnbinom(abundance_intervention_control_sim$total_caught,
               mu=f_control*catch_prop,
               size=phi, log=T))
     intervention_abundance_treatment_likelihood = sum(
       dnbinom(abundance_intervention_treatment_sim$total_caught,
               mu=fc_trt*catch_prop, size=phi, log=T))
     
     # likelihood for blood fed status
     
     baseline_fullyfed_likelihood = sum(dbinom(bloodmeal_baseline_sim$total_full,
                                               size=bloodmeal_baseline_sim$total_evaluated,
                                               prob=prob_inside_fullyfed_control,
                                               log=TRUE))
     
     
     
     baseline_partfed_likelihood = sum(dbinom(bloodmeal_baseline_sim$total_half,
                                              size=bloodmeal_baseline_sim$total_evaluated,
                                              prob=prob_inside_partfed_control,
                                              log=TRUE))
     baseline_unfed_likelihood = sum(dbinom(bloodmeal_baseline_sim$total_empty,
                                            size=bloodmeal_baseline_sim$total_evaluated,
                                            prob=prob_inside_unfed_control,
                                            log=TRUE))
     
     prob_inside_fullyfed_treated[prob_inside_fullyfed_treated == 0] = 0.000000001
     prob_inside_partfed_treated[prob_inside_partfed_treated == 0] = 0.000000001
     prob_inside_unfed_treated[prob_inside_unfed_treated == 0] = 0.000000001
     
     intervention_fullyfed_likelihood = sum(dbinom(bloodmeal_intervention_sim$total_full,
                                                   size=bloodmeal_intervention_sim$total_evaluated,
                                                   prob=prob_inside_fullyfed_treated[bloodmeal_intervention_sim$Cluster],
                                                   log=TRUE))
     intervention_partfed_likelihood = sum(dbinom(bloodmeal_intervention_sim$total_half,
                                                  size=bloodmeal_intervention_sim$total_evaluated,
                                                  prob=prob_inside_partfed_treated[bloodmeal_intervention_sim$Cluster],
                                                  log=TRUE))
     intervention_unfed_likelihood = sum(dbinom(bloodmeal_intervention_sim$total_empty,
                                                size=bloodmeal_intervention_sim$total_evaluated,
                                                prob=prob_inside_unfed_treated[bloodmeal_intervention_sim$Cluster],
                                                log=TRUE))
     
     return(-(epi_likelihood + baseline_parity_likelihood + intervention_parity_likelihood +
                intervention_abundance_control_likelihood + baseline_abundance_likelihood+
                intervention_abundance_treatment_likelihood + 
                baseline_fullyfed_likelihood + baseline_partfed_likelihood +
                baseline_unfed_likelihood + intervention_fullyfed_likelihood + 
                intervention_partfed_likelihood +
                intervention_unfed_likelihood))
     
   }
  
  
  
  
   optimized = mle2(minuslogl = likelihood_sim,
                    start=as.list(c(alpha=alpha, qt=qt, gt=gt, mu=mu, rho=rho, lambda=lambda, af=af, ap=ap, qu=qu,
                                    c=c, gu=gu, n=n, tau=tau, b=b,
                                    X_1=X_1, X_2=X_2, X_3=X_3, X_4=X_4, X_5=X_5, X_6=X_6, X_7=X_7, X_8=X_8, X_9=X_9,
                                    X_10=X_10, X_11=X_11, X_12=X_12,
                                    X_13=X_13, X_14=X_14, X_15=X_15, X_16=X_16, X_17=X_17, X_18=X_18, X_19=X_19,
                                    X_20=X_20, X_21=X_21,
                                    X_22=X_22, X_23=X_23, X_24=X_24, X_25=X_25, X_26=X_26, gtau=gtau, phi=phi,
                                    d=d, catch_prop=catch_prop)), method="L-BFGS-B",
                    upper=upper,lower=lower)
  
  
  
  # return: data frame with a column labeled true and a column labelled inferred
  return(data.frame(true=c(alpha, qt, gt, mu, rho, lambda, af, ap, qu, c, gu, n, tau, b, 
                    X_1, X_2, X_3, X_4, X_5, X_6, X_7, X_8, X_9, X_10, X_11, X_12, 
                    X_13, X_14, X_15, X_16, X_17, X_18, X_19, X_20, X_21, 
                    X_22, X_23, X_24, X_25, X_26, gtau, phi, d, catch_prop),
             inferred=coef(optimized)))
}




#### parallel implementation ####
cl = makeCluster(20)
registerDoParallel(cl)


result = foreach(i=1:nrow(prior.samps), .errorhandling = 'remove',
                 .export=ls(globalenv()),
                 .packages = c("dplyr","mvtnorm","doParallel", "deSolve",
                               "msm", "markovchain", "bbmle")) %dopar% {
                                 rbind(with(prior.samps[i,], find_max_likelihood(X_1, X_2, X_3, X_4, X_5, X_6,
                                                                           X_7, X_8, X_9, X_10,
                                                                           X_11, X_12, X_13, X_14, X_15,
                                                                           X_16, X_17, X_18, X_19, X_20,
                                                                           X_21, X_22, X_23, X_24, X_25,
                                                                           X_26,af,ap,gu ,
                                                                           gtau , qu, 
                                                                           tau ,n,b,
                                                                           c,alpha, mu,
                                                                           gt,rho,
                                                                           qt,lambda,d,
                                                                           catch_prop,phi)), c(1,1))
                                 
                                }
stopCluster(cl)
save(result, file="SimulationResult.RData")