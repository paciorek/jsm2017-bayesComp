### Hyperparameters
sigma2 = 1^2
tau2   = 2^2
mm0    = 0
CC20   = 4


### Generating data
TT     = 150
theta.t = rep(0, TT+1)
y       = rep(0, TT)

theta.t[1] = 2
for(i in 1:TT){
  theta.t[i+1] = rnorm(1, theta.t[i], sqrt(tau2))
  y[i] = rnorm(1, theta.t[i+1], sqrt(sigma2))
}


### Exact solution
mm  = rep(0,TT+1)
CC2 = rep(0,TT+1)
mm[1]  = mm0
CC2[1] = CC20
for(i in 1:TT){
  mm[i+1]  = (y[i]/sigma2 + mm[i]/(tau2 + CC2[i]))/(1/sigma2 + 1/(tau2 + CC2[i]))
  CC2[i+1] = 1/(1/sigma2 + 1/(tau2 + CC2[i]))
}


### SIR filter
BB = 10000
theta.out   = array(0, dim=c(TT+1, BB))
weights.out = array(0, dim=c(TT+1, BB))
ESS         = rep(0, TT)
theta.out[1,]     = rnorm(BB, mm0, sqrt(CC20))
weights.out[1,]   = rep(1/BB, BB)
for(i in 1:TT){
  theta.tilde = sample(theta.out[i,], BB, replace=TRUE, weights.out[i,])
  theta.hat   = rnorm(BB, theta.tilde, sqrt(tau2))
  lweights    = dnorm(y[i], theta.hat, sqrt(sigma2), log=T)
  lweights    = lweights - max(lweights)

  theta.out[i+1,]   = theta.hat
  weights.out[i+1,] = exp(lweights)/sum(exp(lweights))
  ESS[i]            = (sum(weights.out[i+1,])^2)/sum(weights.out[i+1,]^2)
}

mm.SIR = apply(weights.out*theta.out, 1,sum)
CC2.SIR = apply(weights.out*theta.out^2, 1,sum) - mm.SIR^2




## Expected values
quartz()
plot(theta.t, type="n", xlab="t", ylab=expression(paste("E(",theta[t],"|",y[1:t],")")), lwd=2)
lines(mm[-1], col="red", lwd=2, lty=2)
lines(mm.SIR[-1], col="blue", lwd=2, lty=3)
legend(120,4,c("Exact","SIR"), col=c("red","blue"), lty=c(2,3), lwd=2, bty="n")
dev.print(dev=pdf,file="Gaussian_means.pdf")


## Confidence intervals
quartz()
plot(theta.t, type="n", xlab="t", ylab=expression(paste("95% credible interval for ",theta[t])), lwd=2)
lines(mm[-1]-2*sqrt(CC2[-1]), col="red", lwd=1, lty=2)
lines(mm[-1]+2*sqrt(CC2[-1]), col="red", lwd=1, lty=2)
lines(mm.SIR[-1]-2*sqrt(CC2.SIR[-1]), col="blue", lwd=1, lty=3)
lines(mm.SIR[-1]+2*sqrt(CC2.SIR[-1]), col="blue", lwd=1, lty=3)
legend(120,4,c("Exact","SIR"), col=c("red","blue"), lty=c(2,3), lwd=2, bty="n")
dev.print(dev=pdf,file="Gaussian_confint.pdf")


## ESS as a function of time
quartz()
plot(ESS, type="l", ylab="ESS", xlab="t",ylim=c(0,BB))
dev.print(dev=pdf,file="Gaussian_ESSintime.pdf")


## ESS vs logarithm of the variance
quartz()
plot(log(apply(weights.out,1,var)[-1]), ESS, ylab="ESS", xlab=expression(paste("log(Var(",omega[t]^{(b)},"))")),ylim=c(0,BB))
dev.print(dev=pdf,file="Gaussian_ESSvslogVar.pdf")












