rm(list=ls())
library(MCMCpack)
set.seed(81196)  # So that results are reproducible

nests = read.csv("nestsize.csv",header=F)
x = nests$V1

## Initialize the parameters
KK    = 2
w     = 1/2                         
lam   = mean(x)
n     = length(x)

## The actual MCMC algorithm starts here
# Priors
aa  = rep(1,KK)  # Uniform prior on w
eta = mean(x)    # Mean of x data for the prior on mu_k
alpha = 1
beta = 1

# Number of iterations of the sampler
rrr   = 6000
burn  = 1000


# Storing the samples
cc.out    = array(0, dim=c(rrr, n))
w.out     = rep(0, rrr)
lam.out   = rep(0, rrr)
logpost   = rep(0, rrr)

# MCMC iterations
for(s in 1:rrr){
  # Sample the indicators
  cc = rep(0,n)
  for(i in 1:n){
    v = rep(0,KK)
    v[1] = log(w) + ifelse(x[i] == 0, 0,-Inf)  
    v[2] = log(1-w) + dpois(x[i], lam, log=TRUE)
    v = exp(v - max(v))/sum(exp(v - max(v)))
    cc[i] = sample(1:KK, 1, replace=TRUE, prob=v)
  }
  
  # Sample the weights
  w = rbeta(1, aa[1] + sum(cc==1), aa[2] + sum(cc==2))

  # Update Poisson rate
  x_poisson = x[cc == 2]
  lam = rgamma(1, shape = alpha + sum(x_poisson), rate = beta + length(x_poisson))

  # Store samples
  cc.out[s, ]   = cc
  w.out[s]      = w
  lam.out[s] = lam
  
  # Log-posterior calculation
  for(i in 1:n){
    if(cc[i] == 1){
      logpost[s] = logpost[s] + log(w)
    } else {
      logpost[s] = logpost[s] + log(1 - w) + dpois(x[i], lam, log = TRUE)
    }
  }
  
  logpost[s] = logpost[s] + dbeta(w, aa[1], aa[2], log = TRUE)
  logpost[s] = logpost[s] + dgamma(lam, alpha, beta, log = TRUE)

  if(s/500==floor(s/500)){
    print(paste("s =",s))
  }
}

## Plot the logposterior distribution for various samples
par(mfrow=c(1,1))
par(mar=c(4,4,1,1)+0.1)
plot(logpost, type="l", xlab="Iterations", ylab="Log posterior")

xx = 0:max(x)
density.posterior = array(0, dim=c(rrr-burn,length(xx)))
for(s in 1:(rrr - burn)){
  density.posterior[s, ] = w.out[s + burn] * (xx == 0) + 
                           (1 - w.out[s + burn]) * dpois(xx, lam.out[s + burn])
}
density.posterior.m = apply(density.posterior , 2, mean)
density.posterior.lq = apply(density.posterior, 2, quantile, 0.025)
density.posterior.uq = apply(density.posterior, 2, quantile, 0.975)

par(mfrow=c(1,1))
par(mar=c(4,4,1,1)+0.1)
plot(xx, density.posterior.m, type="n",ylim=c(0,max(density.posterior.uq)), xlab="x", ylab="Density")
polygon(c(xx,rev(xx)), c(density.posterior.lq, rev(density.posterior.uq)), col="grey", border="grey")
lines(xx, density.posterior.m, lwd=2)
points(x, rep(0,n), col=cc.out[rrr, ])