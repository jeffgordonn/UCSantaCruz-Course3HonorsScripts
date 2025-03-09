rm(list=ls())
### Loading data and setting up global variables
library(MASS)
library(MCMCpack)
data("faithful")
kk = 2
x = faithful$eruptions
n = length(x)
plot(density(x))
points(x,rep(0,n))

## priors
aa  = rep(1:kk)  
eta = mean(x)    
tau = sqrt(var(x))
dd  = 2
qq  = var(x)/kk

## initialize the parameters
w     = 1/2
mu    = rnorm(kk, mean(x), sd(x))
sigma = rep(sd(x)/kk,kk)
cc    = sample(1:kk, n, replace=T, prob=c(w,1-w))

## set iters
rrr   = 12000
burn  = 2000

## store samples
cc.out    = array(0, dim=c(rrr, n))
w.out     = rep(0, rrr)
mu.out    = array(0, dim=c(rrr, 2))
sigma.out = array(0, dim=c(rrr, 2))
logpost   = rep(0, rrr)

for(s in 1:rrr){
  ## e-step
  for(i in 1:n){
    v    = rep(0,2)
    v[1] = log(w) + dnorm(x[i], mu[1], sigma[1], log=TRUE)  #Compute the log of the weights
    v[2] = log(1-w) + dnorm(x[i], mu[2], sigma[2], log=TRUE)  #Compute the log of the weights
    v    = exp(v - max(v))/sum(exp(v - max(v)))
    cc[i] = sample(1:2, 1, replace=TRUE, prob=v)
  }
  ## m-step
  # weights
  w = rbeta(1, aa[1] + sum(cc==1), aa[2] + sum(cc==2))
  
  for(k in 1:kk){
    # mean
    nk    = sum(cc==k)
    xsumk = sum(x[cc==k])
    tau2.hat = 1/(nk/sigma[k]^2 + 1/tau^2)
    mu.hat  = tau2.hat*(xsumk/sigma[k]^2 + eta/tau^2)
    mu[k]   = rnorm(1, mu.hat, sqrt(tau2.hat))
    # variance
    dd.star = dd + nk/2
    qq.star = qq + sum((x[cc==k] - mu[k])^2)/2
    sigma[k] = sqrt(1/rgamma(1, dd.star, qq.star))
  }
  
  ## samples
  cc.out[s,]    = cc
  w.out[s]     = w
  mu.out[s,]    = mu
  sigma.out[s,] = sigma
  logpost[s] = 0
  for(i in 1:n){
    if(cc[i] == 1){
      logpost[s] = logpost[s] + log(w) + dnorm(x[i], mu[1], sigma[1], log=TRUE)
    }else{
      logpost[s] = logpost[s] + log(1-w) + dnorm(x[i], mu[2], sigma[2], log=TRUE)
    }
  }
  logpost[s] = logpost[s] + dbeta(w, aa[1], aa[2], log=TRUE)
  for(k in 1:2){
    logpost[s] = logpost[s] + dnorm(mu[k], eta, tau, log = T) + dgamma(1/sigma[k]^2, dd, qq)/sigma[k]^4
  }
  if(s/500==floor(s/500)){
    print(paste("s =",s))
  }
}

## Compute the samples of the density over a dense grid
xx  = seq(0,7,length=150)
density.mcmc = array(0, dim=c(rrr-burn,length(xx)))
for(s in 1:(rrr-burn)){
  density.mcmc[s,] = density.mcmc[s,] + w.out[s+burn]*dnorm(xx,mu.out[s+burn,1],sigma.out[s+burn,1]) + 
    (1-w.out[s+burn])*dnorm(xx,mu.out[s+burn,2],sigma.out[s+burn,2])
}
density.mcmc.m = apply(density.mcmc , 2, mean)

## Credible plot
density.mcmc.lq = apply(density.mcmc, 2, quantile, 0.025)
density.mcmc.uq = apply(density.mcmc, 2, quantile, 0.975)
plot(xx, density.mcmc.m, type="n",ylim=c(0,max(density.mcmc.uq)),xlab="Eruptions", ylab="Density")
polygon(c(xx,rev(xx)), c(density.mcmc.lq, rev(density.mcmc.uq)), col="grey", border="grey")
lines(xx, density.mcmc.m, col="black", lwd=2)
points(x, rep(0,n))