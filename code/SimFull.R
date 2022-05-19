

library(tidyverse)
library(fitdistrplus)
library(BayesianTools)
library(postpack)
library(coda)
library(doParallel)

#### load functions and trial data ####

source("functions.R")

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



#### function to perform inference on a data set ####

infer_parameters = function(simulated_data,  index){
  
  # unpack the simulated data
  sc_control_sim = simulated_data[[1]][[1]]
  sc_treatment_sim = simulated_data[[1]][[2]]
  
  parity_baseline_sim = simulated_data[[2]][[1]]
  parity_intervention_control_sim = simulated_data[[2]][[2]]
  parity_intervention_treatment_sim = simulated_data[[2]][[3]]
  
  catch_baseline_sim = simulated_data[[3]][[1]]
  catch_intervention_control_sim = simulated_data[[3]][[2]]
  catch_intervention_treatment_sim = simulated_data[[3]][[3]]
  
  feeding_baseline_sim = simulated_data[[4]][[1]]
  feeding_intervention_control_sim = simulated_data[[4]][[2]]
  feeding_intervention_treatment_sim = simulated_data[[4]][[3]]
  
  likelihood = function(params){
    
    lambda = exp(params[1])
    qu=exp(params[7])
    q_mult = exp(params[5])
    g_mult =exp(params[4])
    gu = sigmoid(params[12])
    au=exp(params[2])
    tau=exp(params[6])
    d = exp(params[9])
    
    rho=sigmoid(params[14])
    alpha=exp(params[3])
    b=sigmoid(params[15])
    c=sigmoid(params[13])
    X= sigmoid(params[11])
    phi=exp(params[10])
    catch_prop = sigmoid(params[16])
    n=exp(params[8])
    
    qt = qu * q_mult
    gt = gu * g_mult
    
    locations_treated = location_probability(au, qu, tau, d, rho, q_mult, alpha, C)
    locations_untreated = location_probability(au, qu, tau, d, rho, q_mult, alpha, C=0)
    
    prop_untreated_treated = sum(locations_treated[1:2])
    prop_transition_treated = sum(locations_treated[3:4])
    prop_treated_treated = sum(locations_treated[5:6])
    
    prop_untreated_untreated = sum(locations_untreated[1:2])
    prop_transition_untreated = sum(locations_untreated[3:4])
    
    prop_unfed_treated = locations_treated[6] / prop_treated_treated
    prop_fed_treated = locations_treated[5] / prop_treated_treated
    
    prop_unfed_untreated = locations_untreated[2] / prop_untreated_untreated
    prop_fed_untreated = locations_untreated[1] / prop_untreated_untreated
    
    prop_unfed_transition_untreated = locations_untreated[4] / prop_transition_untreated
    prop_unfed_untreated_treated = locations_treated[2] / prop_untreated_treated
    prop_unfed_transition_treated = locations_treated[4] / prop_transition_treated
    
    # modified death rate
    gc_treatment = gc_calc(gu, g_mult,
                           pi.tau=prop_transition_treated, pi.U=prop_untreated_treated,
                           pi.T=prop_treated_treated)
    gc_control = gc_calc(gu, g_mult,
                         pi.tau=prop_transition_untreated,
                         pi.U=prop_untreated_untreated,
                         pi.T=0)
    
    # modified biting rates
    ac_control=ac_calc(au, alpha, d, pi.UN=locations_untreated[2], pi.TN=0,
                       pi.tauN = locations_untreated[4])
    ac_treatment = ac_calc(au, alpha, d, pi.UN=locations_treated[2],
                           pi.TN=locations_treated[6],
                           pi.tauN = locations_treated[4])
    
    
    # ratio of mosquitoes to humans
    m_treated = lambda / gc_treatment
    m_control = lambda / gc_control
    
    # expected catch number in houses
    f_treated = m_treated * p * prop_treated_treated
    f_control = m_control * p * prop_untreated_untreated
    
    # parity rates
    parity_rate_control = parity(ac_control, gc_control)
    parity_rate_treatment = parity(ac_treatment, gc_treatment)
    
    # force of infection
    foi_c = force_of_infection(b, lambda, ac_control, c, X, gc_control, n)
    foi_t = force_of_infection(b, lambda, ac_treatment, c, X, gc_treatment, n)
    
    prob_not_sc_control = exp(-foi_c * sc_control_sim$followup_days)
    prob_sc_control =1 - prob_not_sc_control
    
    prob_not_sc_treatment = exp(-foi_t * sc_treatment_sim$followup_days)
    prob_sc_treatment =1 - prob_not_sc_treatment
    
    
    # take care of some zeroes
    
    if(is.nan(prop_unfed_treated)){
      prop_unfed_treated = .0000000001
    }
    
    if(is.nan(prop_fed_treated)){
      prop_fed_treated = .0000000001
    }
    
    if(prop_unfed_treated == 0){
      prop_unfed_treated = .0000000001
    }
    if(prop_fed_treated == 0){
      prop_fed_treated = .0000000001
    }
    
    if(prop_unfed_untreated == 0){
      prop_unfed_untreated = .0000000001
    }
    if(prop_fed_untreated == 0){
      prop_fed_untreated = .0000000001
    }
    
    
    if(f_treated == 0){
      f_treated = .0000000001
    }
    
    if(f_control == 0){
      f_untreated = .0000000001
    }
    
    
    
    prob_sc_control[prob_sc_control == 0] = 0.00000001
    prob_sc_treatment[prob_sc_treatment == 0] = 0.00000001
    
    prob_sc_control[prob_sc_control == 1] = .999999999
    prob_sc_treatment[prob_sc_treatment == 1] = .999999999
    
    bloodmeal_likelihood_intervention_treated = sum(dbinom(feeding_intervention_treatment_sim,
                                                           size=sum(bloodmeal_intervention$total_evaluated[
                                                             bloodmeal_intervention$cluster_trt == "T"]),
                                                           prob=prop_fed_treated, log=T))
    
    bloodmeal_likelihood_intervention_control = sum(dbinom(feeding_intervention_control_sim,
                                                           size=sum(bloodmeal_intervention$total_evaluated[
                                                             bloodmeal_intervention$cluster_trt == "C"]),
                                                           prob=prop_fed_untreated, log=T))
    
    bloodmeal_likelihood_baseline = sum(dbinom(feeding_baseline_sim,
                                               size=sum(bloodmeal_baseline$total_evaluated),
                                               prob=prop_fed_untreated, log=T))
    
    abundance_likelihood_intervention_treated = sum(dnbinom(catch_intervention_treatment_sim,
                                                            mu=f_treated*catch_prop,
                                                            size=phi, log=T))
    abundance_likelihood_intervention_control = sum(dnbinom(catch_intervention_control_sim,
                                                            mu=f_control*catch_prop,
                                                            size=phi, log=T))
    
    abundance_likelihood_baseline = sum(dnbinom(catch_baseline_sim,
                                                mu=f_control*catch_prop,
                                                size=phi, log=T))
    
    parity_likelihood_baseline = sum(dbinom(parity_baseline_sim, 
                                            size=sum(parity_baseline$total_caught),
                                            prob=parity_rate_control, log=T))
    
    parity_likelihood_intervention_control = sum(dbinom(parity_intervention_control_sim, 
                                                        size=sum(parity_intervention$total_caught[parity_intervention$cluster_trt == "C"]),
                                                        prob=parity_rate_control, log=T))
    parity_likelihood_intervention_treatment = sum(dbinom(parity_intervention_treatment_sim, 
                                                          size=sum(parity_intervention$total_caught[parity_intervention$cluster_trt == "T"]),
                                                          prob=parity_rate_treatment, log=T))
    
    epi_likelihood_baseline = sum(dbinom(sc_control_sim$outcome, size=1,
                                         prob=prob_sc_control, log=T))
    epi_likelihood_treatment = sum(dbinom(sc_treatment_sim$outcome, size=1, 
                                          prob=prob_sc_treatment, log=T))
    
    return((bloodmeal_likelihood_baseline + bloodmeal_likelihood_intervention_treated +
              bloodmeal_likelihood_intervention_control +
              abundance_likelihood_intervention_treated + abundance_likelihood_intervention_control +
              abundance_likelihood_baseline +
              parity_likelihood_baseline + parity_likelihood_intervention_control +
              parity_likelihood_intervention_treatment +
              epi_likelihood_baseline + epi_likelihood_treatment))
    
  }
  
  bayesianSetup = createBayesianSetup(likelihood, prior_density, prior_sampler)
  
  convergence = F
  
  output= runMCMC(bayesianSetup, settings=list(iterations = 1000),
                  sampler="DEzs")
  total_iterations_so_far = 1000
  
  while(!convergence){
    
    output= runMCMC(output, settings=list(iterations=10000))
    gelman_output = gelmanDiagnostics(output)
    
    
    if (max(gelman_output$psrf[,2])< 1.1){
      convergence = T
      print("We have converged")
      
      result = getSample(output, start = 1)
      result = tail(result, 1000)
      result = as.data.frame(result)
      names(result) =c("lambda", "au", "a.mult", "g.mult", "q.mult",
                       "tau", "qu", "n","d","phi",
                       "X", "gu", "c", "rho", "b", "catch_prop")
      
      # convert to natural units
      result_exp = exp(result)
      result_sigmoid = sigmoid(result)
      result = cbind(result_exp[,1:10], result_sigmoid[,11:16])
      names(result) =c("lambda", "au", "a.mult", "g.mult", "q.mult",
                       "tau", "qu", "n","d","phi",
                       "X", "gu", "c", "rho","b", "catch_prop")
      
      total_iterations_so_far = total_iterations_so_far + 10000
      
    } else{
      total_iterations_so_far = total_iterations_so_far + 10000
      print(paste("Did not converge on", total_iterations_so_far, "iterations"))
      print(max(gelman_output$psrf[,2]))
      
      
    }
  }
  
  return(list(result, total_iterations_so_far))
}

#### set up and run ####

# import simulated datasets (big file)
load("20SimulatedDatasets.RData")

sigmoid = function(x){
  1/(1+exp(-x))
}


inverse_sigmoid = function(x){
  log(x / (1-x))
}


prior_density = function(params){
  lambda_dens = dunif(exp(params[1]), 0, 100, log=T)
  au_dens = dgamma(exp(params[2]),3, 4, log=T)
  a.mult_dens = dgamma(exp(params[3]), 100, 100*(4/3), log=T)
  g.mult_dens = dgamma(exp(params[4]), 100, 75, log=T)
  q.mult_dens = dgamma(exp(params[5]), 15, 15.4, log=T)
  tau_dens = dgamma(exp(params[6]),  2*3.3, 4, log=T)
  qu_dens = dgamma(exp(params[7]), 2, 4, log=T)
  n_dens = dgamma(exp(params[8]), 40, 40/14, log=T)
  d_dens =  dunif(exp(params[9]), .2, 1, log=T)
  phi_dens = dexp(exp(params[10]),1, log=T)
  X_dens = dbeta(sigmoid(params[11]), 1, 7, log=T)
  gu_dens =  dbeta(sigmoid(params[12]), 4, 18, log=T)
  c_dens = dbeta(sigmoid(params[13]), 1, 1, log=T)
  rho_dens = dbeta(sigmoid(params[14]), 2, 8, log=T)
  b_dens = dbeta(sigmoid(params[15]), 1, 1, log=T)
  catch_prop_dens = dbeta(sigmoid(params[16]), 1, 13/87, log=T)
  
  return(lambda_dens + au_dens + a.mult_dens + g.mult_dens + q.mult_dens + tau_dens+
           n_dens + d_dens + phi_dens + X_dens + gu_dens + qu_dens + c_dens + rho_dens+
           b_dens + catch_prop_dens)
}


prior_sampler = function(n=1){
  lambda_samp = log(runif(n, 0, 100))
  au_samp = log(rgamma(n, 3, 4))
  a.mult_samp = log(rgamma(n, 100, 100*(4/3)))
  g.mult_samp = log(rgamma(n, 100, 75))
  q.mult_samp = log(rgamma(n, 15, 15.4))
  tau_samp = log(rgamma(n, 2*3.3, 4))
  qu_samp = log(rgamma(n, 2, 4))
  n_samp = log(rgamma(n, 40, 40/14))
  d_samp =  log(runif(n, .2, 1))
  phi_samp = log(rexp(n,1))
  X_samp = inverse_sigmoid(rbeta(n, 1, 7))
  gu_samp =  inverse_sigmoid(rbeta(n, 4, 18))
  c_samp  = inverse_sigmoid(rbeta(n, 1, 1))
  rho_samp = inverse_sigmoid(rbeta(n, 2, 8))
  b_samp = inverse_sigmoid(rbeta(n, 1, 1))
  catch_prop_samp = rbeta(n, 1, 13/87)
  if(catch_prop_samp > .99){
    catch_prop_samp = inverse_sigmoid(.99)
  } else{
    catch_prop_samp = inverse_sigmoid(catch_prop_samp)
  }
  return(cbind(lambda_samp, au_samp, a.mult_samp, g.mult_samp, q.mult_samp,
               tau_samp, qu_samp, n_samp, d_samp, phi_samp, X_samp, gu_samp, c_samp,
               rho_samp, b_samp, catch_prop_samp))
}





#### run the algorithm on each simulated data set ####
number_cores = 20
registerDoParallel(cores=number_cores)
cl=makeCluster(20)
result = foreach(i=1:20,.errorhandling = "pass") %dopar% infer_parameters(simulated_data= simulated_datasets[[i]],
                                                                          index=i)
stopCluster(cl)

save(result, file=paste0("SimFull.RData"))