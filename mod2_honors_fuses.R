## creating a 2 component mixture model for EM 
# component 1: exponential distribution
# component 2: log-Gaussian

# linux
# fuses = read.csv("~/Bayesian Statistics/Mixture Models/fuses.csv",header=FALSE)
# windows
rm(list=ls())
fuses = read.csv("./Bayesian Statistics/Mixture Models/fuses.csv",header=FALSE)
set.seed(81196)    # So that results are reproducible (same simulated data every time)

## Run the actual EM algorithm
## Initialize the parameters
w     = 0.1
# exp
rate  = 1/mean(fuses$V1)
# log gaussian
mu    = mean(log(fuses$V1))
sigma = sd(log(fuses$V1))   
# fuses = na.omit(fuses)

s  = 0
sw = FALSE
QQ = -Inf
QQ.out = NULL
epsilon = 10^(-5)


##Checking convergence of the algorithm
while(!sw){
  ## E step
  v = array(0, dim=c(length(fuses$V1),2))
  # exp
  v[,1] = log(w) + dexp(fuses$V1, rate, log=TRUE)
  # log gaussian
  v[,2] = log(1-w) + dlnorm(fuses$V1, mu, sigma, log=TRUE)
  for(i in 1:length(fuses$V1)){
    v[i,] = exp(v[i,] - max(v[i,]))/sum(exp(v[i,] - max(v[i,])))  #Go from logs to actual weights in a numerically stable manner
  }
  
  ## M step
  # Weights 
  w = mean(v[,1])
  # exp
  rate = sum(v[,1]) / sum(v[,1] * fuses$V1)
  # log-norm
  mu = sum(v[,2]*log(fuses$V1))/sum(v[,2])
  sigma = sqrt(sum(v[,2]*(log(fuses$V1) - mu)^2)/sum(v[,2]))
  
  ##Check convergence
  QQn = 0
  for (i in 1:length(fuses$V1)) {
    QQn = QQn + v[i,1] * (log(w) + dexp(fuses$V1[i], rate, log=TRUE)) +
                v[i,2] * (log(1 - w) + dlnorm(fuses$V1[i], mu, sigma, log=TRUE))
  }
  if(abs(QQn-QQ)/abs(QQn)<epsilon){
    sw=TRUE
  }
  QQ = QQn
  QQ.out = c(QQ.out, QQ)
  s = s + 1
  print(paste(s, QQn))
  if(s%%100==0) {
	resp = readline("Continue?: ")
	if (tolower(resp) == "no") {
	  sw=TRUE
	}
  }
}

