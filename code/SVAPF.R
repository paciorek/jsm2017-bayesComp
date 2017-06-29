### Hyperparameters
eta    = log(0.15/sqrt(252))
rho    = 0.95
tau2   = 0.25




### Generating data
TT     = 150
theta.t = rep(0, TT+1)
y       = rep(0, TT)

theta.t[1] = eta
for(i in 1:TT){
  theta.t[i+1] = rnorm(1, eta + rho*(theta.t[i] - eta), sqrt(tau2))
  y[i] = rnorm(1, 0, exp(theta.t[i+1]))
}

quartz()
plot(y, type="l")













### APF filter
BB = 10000
theta.out   = array(0, dim=c(TT+1, BB))
weights.out = array(0, dim=c(TT+1, BB))
ESS.APF     = rep(0, TT)
theta.out[1,]     = rnorm(BB, eta, sqrt(tau2/(1-rho^2)))
weights.out[1,]   = rep(1/BB, BB)
for(i in 1:TT){
  mu          = eta + rho*(theta.out[i,] - eta)
  lvarpi      = log(weights.out[i,]) + dnorm(y[i], 0, exp(mu), log=T)
  varpi       = exp(lvarpi - max(lvarpi))
  varpi       = varpi/sum(varpi)
  ind         = sample(1:BB, BB, replace=TRUE, varpi)
  theta.tilde = theta.out[i,ind]
  theta.hat   = rnorm(BB, eta + rho*(theta.tilde - eta), sqrt(tau2))
  theta.out[i+1,]   = theta.hat
  
  lweights    = dnorm(y[i], 0, exp(theta.hat), log=T) - dnorm(y[i], 0, exp(mu[ind]), log=T)
  lweights    = lweights - max(lweights)
  weights.out[i+1,] = exp(lweights)/sum(exp(lweights))
  ESS.APF[i]        = (sum(weights.out[i+1,])^2)/sum(weights.out[i+1,]^2)
}
mm.APF = apply(weights.out*theta.out, 1,sum)
CC2.APF = apply(weights.out*theta.out^2, 1,sum) - mm.APF^2



### SIR filter
BB = 10000
theta.out   = array(0, dim=c(TT+1, BB))
weights.out = array(0, dim=c(TT+1, BB))
ESS.SIR     = rep(0, TT)
theta.out[1,]     = rnorm(BB, eta, sqrt(tau2/(1-rho^2)))
weights.out[1,]   = rep(1/BB, BB)
for(i in 1:TT){
  varpi       = weights.out[i,]
  ind         = sample(1:BB, BB, replace=TRUE, varpi)
  theta.tilde = theta.out[i,ind]
  theta.hat   = rnorm(BB, eta + rho*(theta.tilde - eta), sqrt(tau2))
  lweights    = dnorm(y[i], 0, exp(theta.hat), log=T)
  lweights    = lweights - max(lweights)
  
  theta.out[i+1,]   = theta.hat
  weights.out[i+1,] = exp(lweights)/sum(exp(lweights))
  ESS.SIR[i]        = (sum(weights.out[i+1,])^2)/sum(weights.out[i+1,]^2)
}
mm.SIR = apply(weights.out*theta.out, 1,sum)
CC2.SIR = apply(weights.out*theta.out^2, 1,sum) - mm.SIR^2






## Expected values
quartz()
plot(theta.t, type="l", xlab="t", ylab=expression(paste("E(",theta[t],"|",y[1:t],")")), lwd=2)
lines(mm.SIR[-1], col="blue", lwd=2, lty=3)
lines(mm.APF[-1], col="green", lwd=2, lty=4)
legend(120,4,c("True","SIR","APF"), col=c("black","blue","green"), lty=c(1,3,4), lwd=2, bty="n")
#dev.print(dev=pdf,file="Gaussian_means_APF.pdf")


## ESS as a function of time
quartz()
plot(ESS.SIR, type="l", ylab="ESS", xlab="t", ylim=c(0,BB), lty=3, col="blue")
lines(ESS.APF,lty=4, col="green")
legend(120,9000,c("SIR","APF"), col=c("blue","green"), lty=c(3,4), lwd=2, bty="n")
#dev.print(dev=pdf,file="Gaussian_ESSintime_APF.pdf")

mean(ESS.SIR)
mean(ESS.APF)












