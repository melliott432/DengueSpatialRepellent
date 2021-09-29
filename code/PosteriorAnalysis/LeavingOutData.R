library(tidyverse)
library(Bolstad2)

#### load all data ####

setwd("~/Dropbox/DengueEntomologicalEffects/results/Posteriors")

load("SequentialMonteCarlo.RData")
samples_AllData = head(samples,1000)

load("SMC_NoAspiration.RData")
samples_NoAspiration = head(samples,1000)

load("SMC_NoParity.RData")
samples_NoParity = head(samples, 1000)

load("SMC_NoSC.RData")
samples_NoSC = head(samples, 1000)

load("SMC_NoBloodfed.RData")
samples_NoBloodfed = head(samples, 1000)


#### Make sure they converge ####

load("SequentialMonteCarlo.RData")
samples_AllData = head(samples,1000)

load("SMC_NoAspiration2.RData")
samples_NoAspiration2 = head(samples,1000)

load("SMC_NoParity2.RData")
samples_NoParity2 = head(samples, 1000)

load("SMC_NoSC2.RData")
samples_NoSC2 = head(samples, 1000)

load("SMC_NoBloodfed2.RData")
samples_NoBloodfed2 = head(samples, 1000)


conv = rep(NA, ncol(samples_NoBloodfed))
for(i in 1:ncol(samples_NoBloodfed)){
  
  theta = cbind(as.matrix(samples_NoBloodfed[,i]), as.matrix(samples_NoBloodfed2[,i]))
  res = GelmanRubin(theta)
  conv[i] = res$R
}

mean(conv > 1.2)

# aspiration hasn't converged
# neither has parity
# neither has SC
# bloodfed is close


#### Load other files ####

setwd("~/Dropbox/DengueEntomologicalEffects/code")
source("functions.R")

#### Mosquito locations ####

mosquito.locations.control.AllData =mosquito.locations.treatment.AllData = 
  mosquito.locations.control.NoAspiration =mosquito.locations.treatment.NoAspiration =
  mosquito.locations.control.NoParity =mosquito.locations.treatment.NoParity =
  mosquito.locations.control.NoSC =mosquito.locations.treatment.NoSC =
  mosquito.locations.control.NoBloodfed =mosquito.locations.treatment.NoBloodfed =
  matrix(rep(NA, nrow(samples_AllData)*9), nrow=nrow(samples_AllData),
                                    ncol=9)

for(i in 1:nrow(samples_AllData)){
  print(i)
  mosquito.locations.control.AllData[i,] = with(samples_AllData[i,],  location_probability(af, ap, qu, tau, d, rho,
                                                                           qt, alpha, C=0))
  mosquito.locations.treatment.AllData[i,] = with(samples_AllData[i,],  location_probability(af, ap,
                                                                             qu, tau,
                                                                             d, rho,
                                                                             qt, alpha,
                                                                             C=.45))
  mosquito.locations.control.NoAspiration[i,] = with(samples_NoAspiration[i,],  location_probability(af, ap, qu, tau, d, rho,
                                                                           qt, alpha, C=0))
  mosquito.locations.treatment.NoAspiration[i,] = with(samples_NoAspiration[i,],  location_probability(af, ap,
                                                                             qu, tau,
                                                                             d, rho,
                                                                             qt, alpha,
                                                                             C=.45))
  
  mosquito.locations.control.NoParity[i,] = with(samples_NoParity[i,],  location_probability(af, ap, qu, tau, d, rho,
                                                                                        qt, alpha, C=0))
  mosquito.locations.treatment.NoParity[i,] = with(samples_NoParity[i,],  location_probability(af, ap,
                                                                                          qu, tau,
                                                                                          d, rho,
                                                                                          qt, alpha,
                                                                                          C=.45))
  
  mosquito.locations.control.NoSC[i,] = with(samples_NoSC[i,],  location_probability(af, ap, qu, tau, d, rho,
                                                                                        qt, alpha, C=0))
  mosquito.locations.treatment.NoSC[i,] = with(samples_NoSC[i,],  location_probability(af, ap,
                                                                                          qu, tau,
                                                                                          d, rho,
                                                                                          qt, alpha,
                                                                                          C=.45))
  
  mosquito.locations.control.NoBloodfed[i,] = with(samples_NoBloodfed[i,],  location_probability(af, ap, qu, tau, d, rho,
                                                                                        qt, alpha, C=0))
  mosquito.locations.treatment.NoBloodfed[i,] = with(samples_NoBloodfed[i,],  location_probability(af, ap,
                                                                                          qu, tau,
                                                                                          d, rho,
                                                                                          qt, alpha,
                                                                                          C=.45))
} 


#### Make a table of estimates ####
number_data = 5
comparison_table = data.frame(alpha=rep(NA, number_data), relative_exit=rep(NA, number_data),
                              mu=rep(NA, number_data),
                              rho=rep(NA, number_data), X_11=rep(NA, number_data),
                              relative_death=rep(NA, number_data), relative_biting=rep(NA, number_data))

df_avail = c("AllData", "NoAspiration", "NoParity", "NoSC", "NoBloodfed")
for(i in 1:length(df_avail)){
  current_samples = get(paste0("samples_",df_avail[i]))
  mosquito.locations.control = get(paste0("mosquito.locations.control.", df_avail[i]))
  mosquito.locations.treatment = get(paste0("mosquito.locations.treatment.", df_avail[i]))
  current_samples$control.pi.U = rowSums(mosquito.locations.control[,1:3])
  current_samples$control.pi.tau = rowSums(mosquito.locations.control[,4:6])
  
  current_samples$treatment.pi.U = rowSums(mosquito.locations.treatment[,1:3])
  current_samples$treatment.pi.tau = rowSums(mosquito.locations.treatment[,4:6])
  current_samples$treatment.pi.T = rowSums(mosquito.locations.treatment[,7:9])
  
  comparison_table[i,] = current_samples %>%
    mutate(ac_control = ac_calc(af, ap, rho, C=0, alpha, qt, qu, tau),
           ac_treatment = ac_calc(af, ap, rho, C=.45, alpha, qt, qu, tau),
           gc_control=gc_calc(mu, qu, qt, gu, gt, gtau, tau, rho, C=0, control.pi.tau, control.pi.U, pi.T=0),
           gc_treatment=gc_calc(mu, qu, qt, gu, gt, gtau, tau, rho, C=.45, treatment.pi.tau, treatment.pi.U, treatment.pi.T)) %>%
    mutate(relative_biting = ac_treatment / ac_control,
           relative_death = gc_treatment/gc_control,
           relative_exit = qt/qu) %>%
    dplyr::select(alpha, relative_exit, mu, rho, X_11, relative_death, relative_biting) %>%
    colMeans() %>%
    as.numeric()
  
}

rownames(comparison_table) = c("All Data", "No Aspiration", "No Parity", "NoSC", "No BloodFed")
colnames(comparison_table) = c("Rel. Biting (house)", "Rel. Exit", "Knockdown", "Repellency", "Prevalence", "Rel. Mortality (Cluster)",
                               "Rel. Biting (Cluster)")
setwd("~/Dropbox/DengueEntomologicalEffects/results/ManuscriptPlots")
write.csv(round(comparison_table, 2), file="LeaveOutDataTypes.csv")     

