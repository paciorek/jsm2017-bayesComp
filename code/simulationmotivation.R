##### Data and basic functions
n1 = 10000
n2 = 100
n3 = 100
x1 = 312
x2 = 9
x3 = 2
predpower = function(theta, alpha, beta){
  eta   = (1-beta)*theta/((1-beta)*theta + alpha*(1-theta))
  return(eta)
}


##### Frequentist solution
# Point estimation
thetah = x1/n1
alphah = x2/n2
betah  = 1 - x3/n3
etah   = predpower(thetah, alphah, betah)
print(paste("MLE for eta is ", etah))

# Interval estimation
B = 10000
x1b = rbinom(B, n1, thetah)
x2b = rbinom(B, n2, alphah)
x3b = rbinom(B, n3, 1-betah)
bootstrapsamples = cbind(x1b/n1, x2b/n2, 1-x3b/n3)
etab = predpower(bootstrapsamples[,1],bootstrapsamples[,2],bootstrapsamples[,3])
print(paste("A 95% confidence interval for eta is (", quantile(etab, 0.025), " , ",quantile(etab, 0.975),")",sep=""))


##### Bayesian solution
# Samples from the posterior distribution
posteriorsamples = cbind(rbeta(B, 0.5+x1, 0.5+n1-x1), rbeta(B, 0.5+x2, 0.5+n2-x2), rbeta(B, 0.5+n3-x3, 0.5+x3))
etap = predpower(posteriorsamples[,1], posteriorsamples[,2], posteriorsamples[,3])
quartz()
hist(etap, main=expression(paste("Posterior distribution for ", eta)), xlab=expression(eta))

# Point estimation
print(paste("Posterior mean for eta ", mean(etap)))
print(paste("Posterior median for eta ", median(etap)))

# Interval estiamtion
print(paste("A 95% confidence interval for eta is (", quantile(etap, 0.025), " , ",quantile(etap, 0.975),")",sep=""))


