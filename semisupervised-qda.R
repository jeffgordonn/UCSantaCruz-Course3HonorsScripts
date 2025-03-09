## Using mixture models for classification in the banknote dataset
## Compare linear and quadratic discriminant analysis and a 
##   (semisupervised) location and scale mixture model with K normals
install.packages("mvtnorm")
library(MASS)
library(mvtnorm)
banknote = load("~/Downloads/banknoteclassification.rdata")
n = dim(banknote.training)[1]  # Size of the training set
m = dim(banknote.test)[1]      # Size of the test set
x = rbind(as.matrix(banknote.training[,-1]), as.matrix(banknote.test[,-1]))   # Create dataset of observations, first n belong to the training set, and the rest belong to the test set
p       = dim(x)[2]              # Number of features
KK      = 3
epsilon = 0.00001


# Initialize the parameters of the algorithm
set.seed(63252)
w   = rep(1,KK)/KK  #Assign equal weight to each component to start with
mu  = rmvnorm(KK, apply(x,2,mean), var(x))   #Cluster centers randomly spread over the support of the data
Sigma      = array(0, dim=c(KK,p,p))  #Initial variances are assumed to be the same
Sigma[1,,] = var(x)/KK  
Sigma[2,,] = var(x)/KK
Sigma[3,,] = var(x)/KK

sw     = FALSE
KL     = -Inf
KL.out = NULL
s      = 0

while (!sw) {
  ## E-step
  v = array(0, dim=c(n+m, KK))
  for (k in 1:KK) {
    v[1:n, k] = ifelse(banknote.training[,1] == k, 0, -Inf)
    v[(n+1):(n+m), k] = log(w[k]) + mvtnorm::dmvnorm(x[(n+1):(n+m), ], mu[k, ], Sigma[k,, ], log=TRUE)
  }
  
  for (i in 1:(n+m)) {
    v[i, ] = exp(v[i, ] - max(v[i, ])) / sum(exp(v[i, ] - max(v[i, ])))
  }
  
  ## M-step
  w = apply(v, 2, mean)
  mu = matrix(0, nrow=KK, ncol=p)
  for (k in 1:KK) {
    mu[k, ] = colSums(v[, k] * x) / sum(v[, k])
  }
  
  Sigma = array(0, dim=c(KK, p, p))
  for (k in 1:KK) {
    for (i in 1:(n+m)) {
      Sigma[k,, ] = Sigma[k,, ] + v[i, k] * (x[i, ] - mu[k, ]) %*% t(x[i, ] - mu[k, ])
    }
    Sigma[k,, ] = Sigma[k,, ] / sum(v[, k]) + diag(1e-6, p)  # Add small regularization
  }
  
  ## Check convergence
  KLn = 0
  for (i in 1:(n+m)) {
    for (k in 1:KK) {
      log_val = log(w[k]) + mvtnorm::dmvnorm(x[i, ], mu[k, ], Sigma[k,, ], log=TRUE)
      if (!is.na(log_val) && !is.nan(log_val) && !is.infinite(log_val)) {
        KLn = KLn + v[i, k] * log_val
      }
    }
  }
  
  if (!is.nan(KLn) && !is.na(KLn) && !is.infinite(KLn) && abs(KLn - KL)/max(abs(KLn), 1e-6) < epsilon) {
    sw = TRUE
  }
  
  KL = KLn
  KL.out = c(KL.out, KL)
  s = s + 1
  print(paste(s, KLn))
}

# Predicted labels
predicted_labels = factor(apply(v, 1, which.max)[(n+1):(n+m)] , levels=1:2, labels=c("genuine", "counterfeit"))
# True labels
true_labels = factor(banknote.test.labels)
# Comparison
sum(predicted_labels != true_labels)

# QDA
modqda = qda(banknote.training[,-1], grouping=banknote.training.labels)
ccpredqda = predict(modqda, newdata=banknote.test[,-1])
sum(ccpredqda$class != true_labels)

print("Semi-Supervised: ")
sum(predicted_labels != true_labels)
print("R QDA")
print(sum(ccpredqda$class != true_labels))
