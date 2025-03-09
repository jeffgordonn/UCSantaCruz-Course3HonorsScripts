library(mvtnorm)

mu1 = c(0,0)
mu2 = c(1/3,1/3)
mu3 = c(-2/3,1/3)
x = c(31,-23)
Sigma = diag(1,2)

l1 = dmvnorm(x,mu1, Sigma, log=T)
l2 = dmvnorm(x,mu2, Sigma, log=T)
l3 = dmvnorm(x,mu3, Sigma, log=T)

exp(l2 - max(l1,l2,l3))/(exp(l1-max(l1,l2,l3)) + exp(l2-max(l1,l2,l3)) + exp(l3-max(l1,l2,l3)))