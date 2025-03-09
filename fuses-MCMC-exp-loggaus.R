# linux
# fuses = read.csv("~/Bayesian Statistics/Mixture Models/fuses.csv",header=FALSE)
# windows
rm(list=ls())
library(MCMCpack)
fuses = read.csv("fuses.csv",header=FALSE)
set.seed(81196)    # So that results are reproducible (same simulated data every time)
x = fuses$V1

## Run the actual EM algorithm
## Initialize the parameters
KK    = 2
n     = length(x)
w     = 0.5
# exp
rate  = 1/mean(fuses$V1)
# log gaussian
mu    = mean(log(fuses$V1))
sigma = sd(log(fuses$V1))  


## The actual MCMC algorithm starts here
# Priors
# weights
aa  = rep(1,KK)
# exponential for IGam
alpha = 2
beta = 1
# log gaussian (gaussian parameters)
eta = 0          
tau = 5
dd  = 2
qq  = 1

# Number of iterations of the sampler
rrr   = 6000
burn  = 1000

# Storing the samples
cc.out    = array(0, dim=c(rrr, n))
w.out     = rep(0, rrr)
mu.out    = rep(0, rrr)
sigma.out = rep(0, rrr)
rate.out   = rep(0, rrr)
logpost   = rep(0, rrr)

# MCMC iterations
for(s in 1:rrr){
  # Sample the indicators
  cc = rep(0,n)
  for(i in 1:n){
    v = rep(0,KK)
    v[1] = log(w) + dexp(x[i], rate, log=TRUE)
    v[2] = log(1-w) +  dlnorm(x[i], mu, sigma, log=TRUE)  
    v = exp(v - max(v))/sum(exp(v - max(v)))
    cc[i] = sample(1:KK, 1, replace=TRUE, prob=v)
  }
  
  # Sample the weights
  w = rbeta(1, aa[1] + sum(cc==1), aa[2] + sum(cc==2))

  ## Exponential
  rate = rgamma(1, shape = alpha + sum(x[cc == 2]), rate = beta + length(x[cc == 2]))

  ## Log Gaussian
  # Sample the means
  for(k in 1:KK){
    nk    = sum(cc==k)
    xsumk = sum(log(x[cc==k]))
    tau2.hat = 1/(nk/sigma^2 + 1/tau^2)
    mu.hat  = tau2.hat*(xsumk/sigma^2 + eta/tau^2)
    mu[k]   = rnorm(1, mu.hat, sqrt(tau2.hat))
  }

  # Sample the variances
  dd.star = dd + n/2
  qq.star = qq + sum((log(x[cc == k]) - mu[k])^2)/2
  sigma = sqrt(rinvgamma(1, dd.star, qq.star))

  # Store samples
  cc.out[s,]   = cc
  w.out[s]     = w
  mu.out[s,]   = mu
  sigma.out[s] = sigma
  rate.out[s]  = rate
  for(i in 1:n){
    if(cc[i]==1){
      logpost[s] = logpost[s] + log(w) + dexp(x[i], rate, log=TRUE)
    }else{
      logpost[s] = logpost[s] + log(1-w) + dlnorm(x[i], mu, sigma, log=TRUE)
    }
  }
  logpost[s] = logpost[s] + dbeta(w, aa[1], aa[2],log = T)
  logpost[s] = logpost[s] + log(dinvgamma(sigma^2, dd, 1/qq))
  logpost[s] = logpost[s] + dnorm(mu, eta, tau, log = T)
  if(s/500==floor(s/500)){
    print(paste("s =",s))
  }
}

## Plot the logposterior distribution for various samples
par(mfrow=c(1,1))
par(mar=c(4,4,1,1)+0.1)
plot(logpost, type="l", xlab="Iterations", ylab="Log posterior")

xx = seq(0,15,length=200)
density.posterior = array(0, dim=c(rrr-burn,length(xx)))
for(s in 1:(rrr-burn)){
  density.posterior[s,] = density.posterior[s,] + w.out[s+burn]*dlnorm(xx,mu.out[s+burn],sigma.out[s+burn]) +
                                                  (1-w.out[s+burn])*dexp(xx,rate.out[s+burn])
}
density.posterior.m = apply(density.posterior , 2, mean)
density.posterior.lq = apply(density.posterior, 2, quantile, 0.025)
density.posterior.uq = apply(density.posterior, 2, quantile, 0.975)
par(mfrow=c(1,1))
par(mar=c(4,4,1,1)+0.1)
plot(xx, density.posterior.m, type="n",ylim=c(0,max(density.posterior.uq)), xlab="x", ylab="Density")
polygon(c(xx,rev(xx)), c(density.posterior.lq, rev(density.posterior.uq)), col="grey", border="grey")
lines(xx, density.posterior.m, lwd=2)
points(x, rep(0,n), col=cc.out[rrr,])