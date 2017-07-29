####### GP example
### Basic functions
rm(list=ls())
library(mvtnorm)



###################################################################################
###################################################################################

### Generate data
d        = 10
MM       = matrix(rnorm(d^2, 0, 1), d, d)
Sigma    = MM%*%t(MM)

### Metropolis-Hastings step
sample.theta = function(Sigma, theta, Omega){
  # This could be improved if we coded the density of the multivariate normal ourselves in terms of Sigma^{-1} and stored it (would save one inversion per iteration)
  d         = length(theta)
  thetap    = as.vector(rmvnorm(1, theta, Omega))
  laceptr   = dmvnorm(thetap, rep(0, d), Sigma, log=T) - dmvnorm(theta, rep(0, d), Sigma, log=T)
  
  
  u = runif(1,0,1)
  if(log(u)<laceptr){
    theta   = thetap
    aceptr   = 1
  }else{
    aceptr   = 0
  }
  return(list(theta=theta, aceptr=aceptr))
}


###################################################################################
###################################################################################

### Run it with a fixed covariance matrix!
repl  = 110000
burn  = 10000

Omega = (2.38^2)*Sigma/d  #Optimal proposal
Omega.d = (0.1^2)*diag(d)/d  #Default proposal

theta = as.vector(rmvnorm(1, rep(0,d), Sigma))


theta.out  = array(0, dim=c(repl-burn, d))
aceptr.mh  = 0

for(s in 1:repl){
  temp   = sample.theta(Sigma, theta, Omega.d)
  theta  = temp$theta
  aceptr.mh = aceptr.mh + temp$aceptr
  if(s/10000==floor(s/10000)){
    print(s)
  }
  if(s > burn){
    theta.out[s-burn,] = theta
  }
}

print(aceptr.mh/repl)
print(paste("Estimated mean of theta_1 is ", mean(theta.out[,1]), "    Estimated variance of theta_1 is ", var(theta.out[,1])))
print(paste("True mean of theta_1 is ", 0, "    True variance of theta_1 is ", Sigma[1,1]))
quartz()
par(mfrow=c(2,2))
hist(theta.out[,1])
plot(theta.out[,1], type="l")
acf(theta.out[,1])
effectiveSize(theta.out[,1])



### Run it now with the adaptive MCMC
repl    = 110000
burn    = 10000
cutoff  = 20000
bet     = 0.05
Omega.d = (0.1^2)*diag(d)/d
theta   = as.vector(rmvnorm(1, rep(0,d), Sigma))


theta.out  = array(0, dim=c(repl-burn, d))
aceptr.mh  = 0

for(s in 1:repl){
  
  if(s<cutoff){
    temp   = sample.theta(Sigma, theta, Omega.d)
  }else{
    u = runif(1)
    if(u < bet){
      temp    = sample.theta(Sigma, theta, Omega.d)
    }else{
      Omega.o = 2.38^2 * var(theta.out[1:(s-1-burn),]) / d
      temp    = sample.theta(Sigma, theta, Omega.o)
    }
  }
  theta  = temp$theta
  aceptr.mh = aceptr.mh + temp$aceptr
  if(s/10000==floor(s/10000)){
    print(s)
  }
  if(s > burn){
    theta.out[s-burn,] = theta
  }
}

print(aceptr.mh/repl)
print(paste("Estimated mean of theta_1 is ", mean(theta.out[,1]), "    Estimated variance of theta_1 is ", var(theta.out[,1])))
print(paste("True mean of theta_1 is ", 0, "    True variance of theta_1 is ", Sigma[1,1]))
quartz()
par(mfrow=c(2,2))
hist(theta.out[,1])
plot(theta.out[,1], type="l")
acf(theta.out[,1])
effectiveSize(theta.out[,1])


###################################################################################
###################################################################################

### The following code generates the graphs included in the slides

###
quartz()
acf(theta.out[,1], main="Default proposal", ylab=expression(paste("ACF of ", theta[1])))
dev.print(dev=pdf, file="adap_acf_def.pdf")
quartz()
plot(theta.out[,1], main="Default proposal", type="l", ylab=expression(paste("Trace of ", theta[1])))
dev.print(dev=pdf, file="adap_trace_def.pdf")
###
quartz()
acf(theta.out[,1], main="Optimal proposal", ylab=expression(paste("ACF of ", theta[1])))
dev.print(dev=pdf, file="adap_acf_opt.pdf")
quartz()
plot(theta.out[,1], main="Optimal proposal", type="l", ylab=expression(paste("Trace of ", theta[1])))
dev.print(dev=pdf, file="adap_trace_opt.pdf")
###
quartz()
acf(theta.out[,1], main="Adaptive proposal", ylab=expression(paste("ACF of ", theta[1])))
dev.print(dev=pdf, file="adap_acf_adp.pdf")
quartz()
plot(theta.out[,1], main="Adaptive proposal", type="l", ylab=expression(paste("Trace of ", theta[1])))
dev.print(dev=pdf, file="adap_trace_adp.pdf")






