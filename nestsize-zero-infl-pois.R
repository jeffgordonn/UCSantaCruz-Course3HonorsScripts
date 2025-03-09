### Zero-inflated poisson
set.seed(81196)    # So that results are reproducible (same simulated data every time)

nests = read.csv("nestsize.csv",header=F)
# head(nests)
# str(nests)

x = nests$V1

# initialization
KK    = 2
n = length(x)
w     = 1/2                         
lam   = mean(x)   
                       
s  = 0
sw = FALSE
QQ = -Inf
QQ.out = NULL
epsilon = 10^(-5)

while(!sw){
  ## E step
  v = array(0, dim=c(n,KK))
  v[,1] = log(w) + ifelse(x == 0, 0,-Inf)
  v[,2] = log(1-w) + dpois(x, lam, log=TRUE)
  for(i in 1:n){
    v[i,] = exp(v[i,] - max(v[i,]))/sum(exp(v[i,] - max(v[i,])))  
  }
  
  ## M step
  # Weights
  w = mean(v[,1])

  # lambda
  lam = sum(v[,2]*x)/sum(v[,2])
  
  # Check convergence
  QQn=sum(v[,1] * log(w) + v[,2] * (log(1-w) + dpois(x, lam, log=TRUE)))
  if(abs(QQn-QQ)/abs(QQn)<epsilon){
    sw=TRUE
  }
  QQ = QQn
  QQ.out = c(QQ.out, QQ)
  s = s + 1
  print(paste(s, QQn))
  
  # Plot current estimate over data
  layout(matrix(c(1,2),2,1), widths=c(1,1), heights=c(1.3,3))
  par(mar=c(3.1,4.1,0.5,0.5))
  plot(QQ.out[1:s],type="l", xlim=c(1,max(10,s)), las=1, ylab="Q", lwd=2)
  
  par(mar=c(5,4,1.5,0.5))
  xx = seq(min(x), max(x), length=200)
  yy = w * (xx == 0) + (1 - w) * dpois(xx, lam)
  hist(x, freq=FALSE, breaks=seq(min(x)-0.5, max(x)+0.5, by=1), main=paste("s =",s,"   Q =", round(QQ.out[s],4)), xlab="x", ylab="Density")
  lines(xx, yy, col="red", lwd=2, lty=2)
  legend("topright",c("Estimate"),col=c("red"), lty=c(2))
}


layout(matrix(c(1,2),2,1), widths=c(1,1), heights=c(1.3,3))
par(mar=c(3.1,4.1,0.5,0.5))
plot(QQ.out[1:s],type="l", xlim=c(1,max(10,s)), las=1, ylab="Q", lwd=2)

par(mar=c(5,4,1.5,0.5))
xx = seq(min(x), max(x), length=200)
yy = w * (xx == 0) + (1 - w) * dpois(xx, lam)
hist(x, freq=FALSE, breaks=seq(min(x)-0.5, max(x)+0.5, by=1), main=paste("s =",s,"   Q =", round(QQ.out[s],4)), xlab="x", ylab="Density")
lines(xx, yy, col="red", lwd=2, lty=2)
legend("topright", c("Estimate"), col=c("red"), lty=c(2), bty="n")

