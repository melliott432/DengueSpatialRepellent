#### Data on parameters 

library(dplyr)
library(fitdistrplus)

# read in data
prior.data = read.csv("PriorData.csv")

# make distributions over each variable
# n a b c r g

# n
n.vals = prior.data %>%
  filter(Parameter == "n") %>%
  pull(Value)

# fit the distribution
n.gamma = fitdist(n.vals, "gamma")$estimate

# a
# mean
quantile(rgamma(1000, 32, 160), c(0.005,0.5, 0.995))
a.gamma = c(16, 160)


# b

b.vals = prior.data %>%
  filter(Parameter == "b") %>%
  pull(Value)
# change 1 to 0.99
b.vals[b.vals == 1] = 0.99
b.beta = fitdist(b.vals, "beta")$estimate

# c

c.vals = prior.data %>%
  filter(Parameter == "c") %>%
  pull(Value)
# change 1 to 0.99
c.vals[c.vals == 1] = 0.99
c.beta = fitdist(c.vals, "beta")$estimate


# r
r.vals = prior.data %>%
  filter(Parameter == "r") %>%
  pull(Value)
r.gamma = fitdist(r.vals, "gamma")$estimate

# g
g.vals = prior.data %>%
  filter(Parameter == "g") %>%
  pull(Value)
g.beta = fitdist(g.vals, "beta")$estimate

# qu and qt 
quantile(rgamma(1000, 12, 40), c(0.005,0.5, 0.995))
q.vals = c(12, 40)

# tau
quantile(rgamma(1000, 16, 32), c(0.005,0.5, 0.995))
tau.vals = c(16, 32)
