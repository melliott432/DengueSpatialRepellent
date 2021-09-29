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
gt = rbeta(n.samps, g.beta[1], g.beta[2])

prior.samples = data.frame(af, ap, b, lambda, c, gu, gt,
                           gtau, n, X)
prior.samples$qu = rgamma(nrow(prior.samples), 12, 40)
prior.samples$tau = rgamma(nrow(prior.samples), 16, 32)

prior.samples$alpha = runif(nrow(prior.samples), 0, 2)
prior.samples$qu = rgamma(nrow(prior.samples), 12, 40)
prior.samples$qt = rgamma(nrow(prior.samples), 12, 40)
prior.samples$gtau = rbeta(nrow(prior.samples), 1.8, 17)
prior.samples$mu = rbeta(nrow(prior.samples), 16, 128)
prior.samples$r = rgamma(nrow(prior.samples), 2.6, 14.6)
prior.samples$tau = rgamma(nrow(prior.samples), 16, 32)
prior.samples$rho = rbeta(nrow(prior.samples), 1.4, 5.6)
prior.samples$d = rgamma(nrow(prior.samples), 32, 160)
prior.samples$phi = runif(nrow(prior.samples), 0, 100)


prior.samples$catch_prop = runif(nrow(prior.samples), 0, 1)
prior.samps = prior.samples

for(i in 1:26){
  prior.samps[,paste0("X_", i)] = runif(nrow(prior.samps), 0, 1)
}

save(prior.samps, file="~/Dropbox/DengueEntomologicalEffects/results/Priors/PriorSamplesRandom.RData")
