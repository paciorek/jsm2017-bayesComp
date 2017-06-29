rm(list=ls())
library(mvtnorm)
library(ellipse)


### Plot of target distribution
Sigma = matrix(c(1,0.95,0.95,1), nrow=2, ncol=2)
mu    = c(15, 45)
e25 = ellipse(Sigma, centre=mu, level = 0.25)
e50 = ellipse(Sigma, centre=mu, level = 0.50)
e80 = ellipse(Sigma, centre=mu, level = 0.80)
e95 = ellipse(Sigma, centre=mu, level = 0.95)
e99 = ellipse(Sigma, centre=mu, level = 0.95)
quartz()
par(mar=c(4,4,1,1)+0.15)
plot(e99, type="n", xlab=expression(theta[1]), ylab=expression(theta[2]))
lines(e25)
lines(e50)
lines(e80)
lines(e95)
dev.print(dev=pdf, file="bivariatenormalhighcor.pdf")


### Hamiltonian method for a bivariate normal
# Test code, bivariate normal, parameters hardcoded
U = function(q)
{
  U     = 0.5*(q - mu)%*%solve(Sigma)%*%t(q - mu)
  #U     = dmvnorm(q, mu, Sigma,log=T)
  return(U)
}


grad_U = function(q)
{
  Sigma  = matrix(c(1,0.95,0.95,1), nrow=2, ncol=2)
  mu     = c(15, 45)
  grad_U = solve(Sigma)%*%t(q - mu)
  return(as.vector(grad_U))
}


HMC = function (U, grad_U, epsilon, L, current_q)
{
  q = current_q
  p = rnorm(length(q),0,1)  # independent standard normal variates
  current_p = p
  # Make a half step for momentum at the beginning
  p = p - epsilon * grad_U(q) / 2
  # Alternate full steps for position and momentum
  for (i in 1:L){
    # Make a full step for the position
    q = q + epsilon * p
    # Make a full step for the momentum, except at end of trajectory
    if (i!=L) p = p - epsilon * grad_U(q)
  }
  # Make a half step for momentum at the end.
  p = p - epsilon * grad_U(q) / 2
  # Negate momentum at end of trajectory to make the proposal symmetric
  p = -p
  # Evaluate potential and kinetic energies at start and end of trajectory
  current_U = U(current_q)
  current_K = sum(current_p^2) / 2
  proposed_U = U(q)
  proposed_K = sum(p^2) / 2
  # Accept or reject the state at end of trajectory, returning either
  # the position at the end of the trajectory or the initial position
  if (runif(1) < exp(current_U-proposed_U+current_K-proposed_K)){
    aceptr = 1
    return(list(aceptr=aceptr, q=q))  # accept
  }else{
    aceptr = 0
    return(list(aceptr=aceptr, q=current_q))  # reject
  }
}

rrr       = 10000
sampler1 = function(rrr, epsilon, L){
  theta.out = array(0, dim=c(rrr,2))
  aceptrt = 0

  theta = rmvnorm(1, c(0, 0), matrix(c(1,0,0,1), nrow=2, ncol=2))
  for(s in 1:rrr){
    temp          = HMC(U, grad_U, epsilon, L, current_q=theta)
    theta         = temp$q
    aceptrt       = aceptrt + temp$aceptr
    theta.out[s,] = theta
    print(s)
  }
  print(paste("Acceptance probability = ", aceptrt/rrr))
  return(theta.out)
}

res1 = sampler1(rrr, epsilon=0.25, L=25)
quartz()
par(mfrow=c(2,1))
par(mar=c(4,4,1,1)+0.15)
plot(res1[,1], ylab=expression(theta[1]), xlab="Iteration", type="l")
abline(h=mu[1], col="grey")
par(mar=c(4,4,1,1)+0.15)
acf(res1[,1], main="")
dev.print(dev=pdf, file="bivariatenormal_hamiltonian.pdf")




##############################################
# Optional Gaussian Random walk Metropolis-Hastings algorithm

sample.theta = function(thetac, Omega){
  thetap = rmvnorm(1, thetac, Omega)
  laceptr = dmvnorm(thetap, mu, Sigma, log=T) - dmvnorm(thetac, mu, Sigma, log=T)
  if (log(runif(1)) < laceptr){
    aceptr = 1
    return(list(aceptr=aceptr, theta=thetap))  # accept
  }else{
    aceptr = 0
    return(list(aceptr=aceptr, theta=thetac))  # reject
  }
}



sampler2 = function(rrr){
  Omega = ((2.38^2)/2)*Sigma
  theta.out = array(0, dim=c(rrr,2))
  aceptrt = 0
  theta = rmvnorm(1, c(0, 0), matrix(c(1,0,0,1), nrow=2, ncol=2))
  for(s in 1:rrr)
  {
    temp          = sample.theta(thetac=theta, Omega)
    theta         = temp$theta
    aceptrt       = aceptrt + temp$aceptr
    theta.out[s,] = theta
    print(s)
  }
  print(paste("Acceptance probability = ", aceptrt/rrr))
  return(theta.out)
}
res2 = sampler2(rrr)


quartz()
par(mfrow=c(2,1))
par(mar=c(4,4,1,1)+0.15)
plot(res2[,1], ylab=expression(theta[1]), xlab="Iteration", type="l")
abline(h=mu[1], col="grey")
par(mar=c(4,4,1,1)+0.15)
acf(res2[,1], main="")
dev.print(dev=pdf, file="bivariatenormal_RW.pdf")



