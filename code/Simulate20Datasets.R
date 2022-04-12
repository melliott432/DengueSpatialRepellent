library(dplyr)
library(fitdistrplus)

#### import trial data ####
setwd("~/Dropbox/DengueEntomologicalEffects/Analyses/BayesianTools/code")
source("functions.R")

setwd("~/Dropbox/DengueEntomologicalEffects/Analyses/BayesianTools/data")
# load trial data 
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


#### function to generate a single data set ####
generate_data = function(au_true, d_true, gu_true, qu_true, tau_true, n_true, b_true, c_true,
                         X_true, a.mult_true, g.mult_true, rho_true, q.mult_true, lambda_true,
                         phi_true, catch_prop_true){
  
  # modified exit and mortality rates
  qt_true = qu_true * q.mult_true
  gt_true = gu_true * g.mult_true
  
  # mosquito locations, control and treatment
  locations_treated_true = location_probability(au_true, qu_true, tau_true, d_true, rho_true,
                                                q.mult_true, a.mult_true, C)
  locations_untreated_true = location_probability(au_true, qu_true, tau_true, d_true, rho_true,
                                                  q.mult_true, a.mult_true, C=0)
  
  # unpack these 
  U_untreated_true = sum(locations_untreated_true[1:2])
  tau_untreated_true = sum(locations_untreated_true[3:4])
  
  U_treated_true = sum(locations_treated_true[1:2])
  tau_treated_true = sum(locations_treated_true[3:4])
  T_treated_true = sum(locations_treated_true[5:6])
  
  # find the proportion of mosquitoes in control homes and fed, and in treatment homes and fed
  prop_fed_control_true = locations_untreated_true[1] / U_untreated_true
  prop_fed_treated_true = locations_treated_true[5] / T_treated_true
  
  # find the modified death rates for control and treated
  overall_death_rate_control_true = gc_calc(gu_true, g.mult_true,
                                            pi.tau=tau_untreated_true,
                                            pi.U=U_untreated_true,
                                            pi.T=0)
  overall_death_rate_treatment_true = gc_calc(gu_true, g.mult_true,
                                              pi.tau=tau_treated_true,
                                              pi.U=U_treated_true,
                                              pi.T=T_treated_true)
  
  # find the modified blood feeding rate for control and treatment
  # why higher biting in treated clusters?
  overall_biting_control_true =ac_calc(au_true, a.mult_true, d_true,
                                       pi.UN=locations_untreated_true[2], pi.TN=0,
                                       pi.tauN = locations_untreated_true[4])
  overall_biting_treated_true = ac_calc(au_true, a.mult_true, d_true,
                                        pi.UN=locations_treated_true[2],
                                        pi.TN=locations_treated_true[6],
                                        pi.tauN = locations_treated_true[4])
  
  # find the force of infection (daily)
  foi_control_true = force_of_infection(b_true, lambda_true, 
                                        a=overall_biting_control_true, c=c_true,
                                        X_true,
                                        g=overall_death_rate_control_true, n_true)
  foi_treatment_true = force_of_infection(b_true, lambda_true,
                                          a=overall_biting_treated_true,
                                          c=c_true,
                                          X_true,
                                          g=overall_death_rate_treatment_true, n_true)
  
  # find the parity rates
  parity_control_true = parity(overall_biting_control_true, overall_death_rate_control_true)
  parity_treatment_true = parity(overall_biting_treated_true, overall_death_rate_treatment_true)
  
  # mosquito to human ratio
  m_control_true = lambda_true / overall_death_rate_control_true
  m_treated_true = lambda_true / overall_death_rate_treatment_true
  
  # find the expected number of mosquitoes located inside a single house
  f_control_true = m_control_true * p * U_untreated_true
  f_treatment_true = m_treated_true * p * T_treated_true
  
  # sample size for epidemiology- vector of followup times
  followup_days_control = epi %>%
    filter(Cluster_Allocation == "C") %>%
    pull(followup_time) *365
  
  followup_days_treatment = epi %>%
    filter(Cluster_Allocation == "T") %>%
    pull(followup_time) *365
  
  # get the probability of seroconverting for each person during this time
  prob_seroconverting_control = 1 - exp(-foi_control_true*followup_days_control)
  prob_seroconverting_treatment = 1 - exp(-foi_treatment_true*followup_days_treatment)
  
  # for each person, get a binomial yes/no outcome
  sc_data_control = data.frame(outcome = rbinom(length(followup_days_control), size=1,
                                                prob=prob_seroconverting_control),
                               followup_days = followup_days_control)
  sc_data_treatment = data.frame(outcome= rbinom(length(followup_days_treatment), size=1,
                                        prob=prob_seroconverting_treatment),
                                 followup_days = followup_days_treatment)
  
  #### generate parity data
  parity_data_baseline = rbinom(n=1, size=sum(parity_baseline$total_caught), prob=parity_control_true)
  parity_data_intervention_control = rbinom(n=1,
                                            size=sum(parity_intervention[parity_intervention$cluster_trt == "C",]$total_caught),
                                            prob=parity_control_true)
  parity_data_intervention_treatment = rbinom(n=1,
                                              size=sum(parity_intervention[parity_intervention$cluster_trt == "T",]$total_caught),
                                              prob=parity_treatment_true)
  
  #### generate Prokopack data
  catch_data_baseline = rnbinom(nrow(abundance_baseline), mu = f_control_true * catch_prop_true, 
                                size=phi_true)
  catch_data_intervention_control= rnbinom(nrow(abundance_intervention_control),
                                           mu = f_control_true * catch_prop_true, 
                                           size=phi_true)
  catch_data_intervention_treatment= rnbinom(nrow(abundance_intervention_treatment),
                                             mu = f_treatment_true * catch_prop_true, 
                                             size=phi_true)
  
  #### generate blood feeding data
  
  feeding_data_baseline = rbinom(1, size=sum(bloodmeal_baseline$total_evaluated),
                                 prob=prop_fed_control_true)
  feeding_data_intervention_control = rbinom(1,
                                             size=sum(bloodmeal_intervention[bloodmeal_intervention$cluster_trt ==
                                                                               "C",]$total_evaluated),
                                             prob=prop_fed_control_true)
  feeding_data_intervention_treatment = rbinom(1,
                                               size=sum(bloodmeal_intervention[bloodmeal_intervention$cluster_trt ==
                                                                                 "T",]$total_evaluated),
                                               prob=prop_fed_treated_true)
  
  # return a list of all data
  return(list(
    # epidemiology 
    list(sc_data_control, sc_data_treatment),
    list(parity_data_baseline, parity_data_intervention_control, parity_data_intervention_treatment),
    list(catch_data_baseline, catch_data_intervention_control, catch_data_intervention_treatment),
    list(feeding_data_baseline, feeding_data_intervention_control, 
         feeding_data_intervention_treatment)
  ))
  
}

#### generate parameter values to generate from ####

N=20000

# These come from ten bosch et al. She only gave default values, so I used these
# as the mean with a pretty narrow 95% CI
au_true_candidates = rgamma(N, 3, 4)
gu_true_candidates = rbeta(N, 4, 18)
qu_true_candidates = rgamma(N, 2, 4)
tau_true_candidates = rgamma(N, 2*3.3, 4)
n_true_candidates = rgamma(N, 40, 40/14)
b_true_candidates = rbeta(N, 1, 1)
c_true_candidates =rbeta(N, 1, 1)
X_true_candidates = rbeta(N, 1, 7)


# these are based on full posterior HDI
a.mult_true_candidates = rgamma(N, 100, 100*(4/3))
g.mult_true_candidates = rgamma(N, 100, 75)
rho_true_candidates = rbeta(N, 2, 8)
q.mult_true_candidates = rgamma(N, 15, 15.4)

# This is informed from "a new cost effective" article on Prokopack
catch_prop_true_candidates = rbeta(N, 1, 13/87)

# These are pretty flat because i don't know
lambda_true_candidates= runif(N, 0, 100)
phi_true_candidates = runif(N, 0, 1)
d_true_candidates = runif(N, .2, 1)


candidate_df = data.frame(au_true=au_true_candidates, d_true=d_true_candidates,
                          gu_true=gu_true_candidates,
                          qu_true=qu_true_candidates, tau_true=tau_true_candidates,
                          n_true=n_true_candidates,
                          b_true=b_true_candidates, 
                          X_true=X_true_candidates, 
                          a.mult_true=a.mult_true_candidates, 
                          g.mult_true=g.mult_true_candidates,
                          rho_true=rho_true_candidates, q.mult_true=q.mult_true_candidates,
                          lambda_true = lambda_true_candidates, c_true = c_true_candidates,
                          phi_true = phi_true_candidates,
                          catch_prop_true = catch_prop_true_candidates)

# for each parameter set, calculate the mosquito locations and overall death and biting rates
candidate_df$foi = 0
for(i in 1:nrow(candidate_df)){
  print(i)
  moz_locations = with(candidate_df[i,], location_probability(au_true, qu_true, tau_true,
                                                              d_true, rho_true,
                                                q.mult_true, a.mult_true, C=0))
  
  overall_death = with(candidate_df[i,], gc_calc(gu_true, g.mult_true,
                                            pi.tau=sum(moz_locations[3:4]),
                                            pi.U=sum(moz_locations[1:2]),
                                            pi.T=0))
  
  overall_biting = with(candidate_df[i,], ac_calc(au_true, a.mult_true, d_true,
                                                  pi.UN=moz_locations[2], pi.TN=0,
                                                  pi.tauN = moz_locations[4]))
  candidate_df$foi[i] = with(candidate_df[i,], force_of_infection(b_true, lambda_true,
                                                                  overall_biting, c_true, X_true,
                                                                  overall_death, n_true))
}


# take only those between 0.2 and 2
combos_to_try = candidate_df[candidate_df$foi > 0 & candidate_df$foi < 0.33/365,]

# run again to generate 1000 more data sets
combos_to_try = tail(combos_to_try, 20)
#### generate 1000 datasets ####



simulated_datasets = list()
for(i in 1:20){
  print(i)
  simulated_datasets[[i]] = generate_data(au_true=combos_to_try$au_true[i],
                                          d_true=combos_to_try$d_true[i],
                                          gu_true=combos_to_try$gu_true[i],
                                          qu_true=combos_to_try$qu_true[i],
                                          tau_true=combos_to_try$tau_true[i],
                                          n_true=combos_to_try$n_true[i],
                                          b_true=combos_to_try$b_true[i],
                                          X_true=combos_to_try$X_true[i], 
                                         a.mult_true=combos_to_try$a.mult_true[i],
                                          g.mult_true=combos_to_try$g.mult_true[i],
                                         rho_true=combos_to_try$rho_true[i],
                                          q.mult_true=combos_to_try$q.mult_true[i],
                                         lambda_true=combos_to_try$lambda_true[i],
                                         c_true = combos_to_try$c_true[i],
                                         phi_true=combos_to_try$phi_true[i],
                                         catch_prop_true = combos_to_try$catch_prop_true[i])
}

save(simulated_datasets, combos_to_try, file="20SimulatedDatasets.RData")


#### marginals ####
N=100
qu_true_candidates = rgamma(N, 2, 4)
tau_true_candidates = rgamma(N, 2*3.3, 4)
n_true_candidates = rgamma(N, 40, 40/14)
b_true_candidates = rbeta(N, 1, 1)
c_true_candidates =rbeta(N, 1, 1)
rho_true_candidates = rbeta(N, 2, 8)
catch_prop_true_candidates = rbeta(N, 1, 13/87)
phi_true_candidates = runif(N, 0, 1)
d_true_candidates = runif(N, .2, 1)



marginal_df = data.frame(d_true=d_true_candidates,
                          qu_true=qu_true_candidates,
                          b_true=b_true_candidates, 
                          rho_true=rho_true_candidates, 
                          phi_true = phi_true_candidates)
