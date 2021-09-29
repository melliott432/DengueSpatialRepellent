#### Load important files ####


source("functions.R")
source("PriorCleaning.R")

#### Generate a prior data set ####
n.samps = 1000
af = rgamma(n.samps, 16, 160)
ap = rgamma(n.samps, 16, 160)
b = rbeta(n.samps, b.beta[1], b.beta[2])
lambda=runif(n.samps, 0, 100)
c = rbeta(n.samps, c.beta[1], c.beta[2])
gu = rbeta(n.samps, g.beta[1], g.beta[2])
n = rgamma(n.samps, n.gamma[1], n.gamma[2])
X = runif(n.samps, 0, .1)
gtau = rbeta(n.samps, g.beta[1], g.beta[2])

# make into a data frame
prior.samples = data.frame(af, ap, b, lambda, c, gu,
                           gtau, n, X)
prior.samples$qu = rgamma(nrow(prior.samples), 12, 40)
prior.samples$tau = rgamma(nrow(prior.samples), 16, 32)

prior.samples$ac_control = with(prior.samples, 
                                ac_calc(af, ap, rho=0, C=0, alpha=1, qt=qu, qu, tau))

prior.samples$d = rgamma(nrow(prior.samples), 32, 160)


location_matrix = matrix(rep(NA, 6*nrow(prior.samples)), 
                         nrow=nrow(prior.samples), ncol=6)
for(i in 1:nrow(prior.samples)){
  location_matrix[i,] = with(prior.samples[i,], location_probability(af, ap, qu, 
                                                                 tau, d, 
                                                                 rho=0, qt=qu,
                                                                 alpha=1, C=0)[1:6])
}

prior.samples$pi.U = rowSums(location_matrix[,1:3])
prior.samples$pi.tau = rowSums(location_matrix[,4:6])

prior.samples$gc_control = with(prior.samples,
                                gc_calc(mu=0, qu, qt=1, gu, gt=1, gtau, tau, rho=0, 
                                        C=0, pi.tau, pi.U, pi.T=0))


# calculate the force of infection in each row
prior.samples$foi = NA

for(i in 1:nrow(prior.samples)){
  prior.samples$foi[i] = with(prior.samples[i,], 
                              force_of_infection(b, lambda, ac_control,
                                                 c, X, gc_control, n))
}

# find minimum and maximum prevalence that also work
prior.samples$minX = NA
prior.samples$maxX = NA

# generate columns that hold prevalences for all 26 clusters
for(i in 1:26){
  prior.samples[paste0("X_", i)] = NA
}

# loop through and mind the minimum and maximum X that works

min_FOI = 0.06
max_FOI = 0.6

# ~40% of samples are kept each time
filter.conditions ={(prior.samples$foi > min_FOI) & (prior.samples$foi < max_FOI)}

prior.samples = prior.samples[filter.conditions,]

#### min max x finder ####
min_max_X_finder = function(df){
  for(i in 1:nrow(df)){
    print(i)

    current.sample = df[i,]
    current.sample["g"] = current.sample["gc_control"]
    
    
    # ONCE for untreated
    min.X = current.sample$X
    max.X = current.sample$X
    reached_min_X = F
    reached_max_X = F
    
    while(reached_min_X == F){
      
      proposed_min_X = min.X - .0001
      
      proposed_min_foi = with(current.sample, force_of_infection(b, lambda, ac_control,
                                                                 c, proposed_min_X,
                                                                 g, n))
      if(proposed_min_foi < min_FOI){
        reached_min_X = T
      } else{
        min.X = proposed_min_X
      }
    }
    df$minX[i] = min.X
    
    max.X = current.sample$X
    reached_max_X = F
    
    while(reached_max_X == F){
      
      proposed_max_X = max.X + .0001
      
      if(proposed_max_X > 1){
        df$maxX[i] = 1
        break
      }
      
      proposed_max_foi = with(current.sample, force_of_infection(b,
                                                                 lambda,
                                                                 ac_control,
                                                                 c, 
                                                                 proposed_max_X,
                                                                 g, n))
      if(proposed_max_foi > max_FOI){
        reached_max_X = T
      } else{
        max.X = proposed_max_X
      }
    }
    df$maxX[i] = max.X
  }
  return(df)
}

#### Apply ####

tester=min_max_X_finder(prior.samples)


cluster_X_generator = function(df){
  for(i in 1:nrow(df)){
    
    current.sample=df[i,]
    
    df[i, 20:45] = runif(26, min=as.numeric(current.sample["minX"]), max=as.numeric(current.sample["maxX"]))
  }
  
  return(df)
}

tester2 = cluster_X_generator(tester)

#### automated procedure ####
prior.samps = tester2

while(nrow(prior.samps) < 100000){
  
  n.samps = 100000
  af = rgamma(n.samps, 16, 160)
  ap = rgamma(n.samps, 16, 160)
  b = rbeta(n.samps, b.beta[1], b.beta[2])
  lambda=runif(n.samps, 0, 100)
  c = rbeta(n.samps, c.beta[1], c.beta[2])
  gu = rbeta(n.samps, g.beta[1], g.beta[2])
  n = rgamma(n.samps, n.gamma[1], n.gamma[2])
  X = runif(n.samps, 0, .1)
  gtau = rbeta(n.samps, g.beta[1], g.beta[2])
  
  # make into a data frame
  prior.samples = data.frame(af, ap, b, lambda, c, gu,
                             gtau, n, X)
  prior.samples$qu = rgamma(nrow(prior.samples), 12, 40)
  prior.samples$tau = rgamma(nrow(prior.samples), 16, 32)
  
  prior.samples$ac_control = with(prior.samples, 
                                  ac_calc(af, ap, rho=0, C=0, alpha=1, qt=qu, qu, tau))
  
  prior.samples$d = rgamma(nrow(prior.samples), 32, 160)
  
  
  location_matrix = matrix(rep(NA, 6*nrow(prior.samples)), 
                           nrow=nrow(prior.samples), ncol=6)
  for(i in 1:nrow(prior.samples)){
    location_matrix[i,] = with(prior.samples[i,], location_probability(af, ap, qu, 
                                                                       tau, d, 
                                                                       rho=0, qt=qu,
                                                                       alpha=1, C=0)[1:6])
  }
  
  prior.samples$pi.U = rowSums(location_matrix[,1:3])
  prior.samples$pi.tau = rowSums(location_matrix[,4:6])
  
  prior.samples$gc_control = with(prior.samples,
                                  gc_calc(mu=0, qu, qt=1, gu, gt=1, gtau, tau, rho=0, 
                                          C=0, pi.tau, pi.U, pi.T=0))
  
  
  # calculate the force of infection in each row
  prior.samples$foi = NA
  
  for(i in 1:nrow(prior.samples)){
    prior.samples$foi[i] = with(prior.samples[i,], 
                                force_of_infection(b, lambda, ac_control,
                                                   c, X, gc_control, n))
  }
  
  # find minimum and maximum prevalence that also work
  prior.samples$minX = NA
  prior.samples$maxX = NA
  
  # generate columns that hold prevalences for all 26 clusters
  for(i in 1:26){
    prior.samples[paste0("X_", i)] = NA
  }
  
  # loop through and mind the minimum and maximum X that works
  
  min_FOI = 0.06
  max_FOI = 0.6
  
  # ~40% of samples are kept each time
  filter.conditions ={(prior.samples$foi > min_FOI) & (prior.samples$foi < max_FOI)}
  
  prior.samples = prior.samples[filter.conditions,]
  
  # find min and max X
  prior.samples = min_max_X_finder(prior.samples)
  
  # generate new values
  prior.samples = cluster_X_generator(prior.samples)
  
  # add the values to prior.samps
  prior.samps = rbind(prior.samps, prior.samples)
  
  
  
  save(prior.samps, file="PriorSamples.RData")
  
}


#### add in other variables ####

prior.samps$alpha = runif(nrow(prior.samps), 0, 2)
prior.samps$qt = rgamma(nrow(prior.samps), 12, 40)
prior.samps$mu = rbeta(nrow(prior.samps), 16, 128)
prior.samps$rho = rbeta(nrow(prior.samps), 1.4, 5.6)
prior.samps$phi = runif(nrow(prior.samps), 0, 100)
prior.samps$catch_prop = runif(nrow(prior.samps), 0, 1)
prior.samps$gt =  rbeta(nrow(prior.samps), g.beta[1], g.beta[2])






save(prior.samps, file="PriorSamples.RData")
