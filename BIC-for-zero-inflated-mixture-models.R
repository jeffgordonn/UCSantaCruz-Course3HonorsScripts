## Full Bayesian estimation of a mixture model for density estimation in the galaxies xset
rm(list=ls())

### Loading x and setting up global variables
nests = read.csv("nestsize.csv",header=T)
x  = nests$X4

poisEM = function(x,max_iter=1000,tol=1e-8) {
  lam = mean(x)
  log_like = sum(dpois(x, lam, log=TRUE))

  twicenegloglik = 2 * (-log_like)
  return(list(lam=lam,twicenegloglik=twicenegloglik))
}

zipEM = function(x,max_iter=1000,tol=1e-8) {
  n = length(x)
  w = 0.1
  lam = mean(x)
  olog_like = 0
  for (iter in 1:max_iter) {
    # e-step
    tau = w* (x==0)/(w*(x==0) + (1-w) * dpois(x,lam))
    # m-step
    w = sum(tau)/n
    lam = sum((1 - tau) * x) / sum(1 - tau)
    log_like =  sum(log(w * (x == 0) + (1 - w) * dpois(x, lam))) 
    if (abs(log_like-olog_like) < tol) {
      break
    }
    olog_like = log_like
  }
  twicenegloglik = 2 * (-log_like)
  return(list(w=w,lam=lam,twicenegloglik=twicenegloglik))
}
### Comparison
n =300
# pois
res1 = poisEM(x)


lam1 = res1$lam
twneglogl1 = res1$twicenegloglik
BIC1 = twneglogl1 + log(n)


# zip
res2 = zipEM(x)

w2 = res2$w
lam2 = res2$lam
twneglogl2 = res2$twicenegloglik
BIC2 = twneglogl2 + log(n)

print(BIC1)
print(BIC2)