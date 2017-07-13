library(nimble)
### Hyperparameters
sigma2 = 1^2
tau2   = 2^2
mm0    = 0
CC20   = 4

### Nimble Gausian SSM Model Code
gaussianSSM.code <- nimbleCode({
  theta.t[1] ~ dnorm(mm0, sd = sqrt(CC20))
  for(i in 1:TT){
    theta.t[i+1] ~ dnorm(theta.t[i], sd = sqrt(tau2))
    y[i] ~ dnorm(theta.t[i+1], sd = sqrt(sigma2))
  }
})


### Generating data
set.seed(1)
TT     = 150
theta.t = rep(0, TT+1)
y       = rep(0, TT)
theta.t[1] = 2
for(i in 1:TT){
  theta.t[i+1] = rnorm(1, theta.t[i], sqrt(tau2))
  y[i] = rnorm(1, theta.t[i+1], sqrt(sigma2))
}

### Creating Nimble Gausian SSM 
gaussianSSM <- nimbleModel(code = gaussianSSM.code, data = list(y = y), 
                           constants = list(sigma2 = sigma2, tau2 = tau2,
                                            mm0 = mm0, CC20 = CC20, TT = TT),
                           dimensions = list(theta.t = TT+1, y = TT))
CgaussianSSM <- compileNimble(gaussianSSM)

### Creating Filter.
## "thresh = 1" option ensures that resampling will happen at each time point.
gaussianSSM.SIRfilter <- buildBootstrapFilter(gaussianSSM, nodes = 'theta.t',
                                 control = list(thresh = 1, saveAll = TRUE))
CgaussianSSM.SIRfilter <- compileNimble(gaussianSSM.SIRfilter,
                                        project = gaussianSSM)
BB = 10000
CgaussianSSM.SIRfilter$run(BB)

### The first 151 columns of the output will be particle samples for all
### 151 time points.
particle.samples <- as.matrix(CgaussianSSM.SIRfilter$mvWSamples)[, 1:151]

### Columns 152 through 302 of the output will be the unnormalized log 
### weights for each particle at each time point.
log.weights.out <- as.matrix(CgaussianSSM.SIRfilter$mvWSamples)[, 152:302]
### Below, we un-log and normalize the weights
weights.out <-  exp(log.weights.out)/apply(exp(log.weights.out), 2, sum)

mm.SIR <- apply(weights.out*particle.samples, 2, sum)
CC2.SIR <- apply(weights.out*particle.samples^2, 2, sum) - mm.SIR^2
ESS <- CgaussianSSM.SIRfilter$returnESS()

### ESS averaged over all time points. Time point 1 corresponds to the 
### prior distribution placed on theta in the model, so we drop it
print(mean(ESS[-1]))

### Exact solution
mm  = rep(0,TT+1)
CC2 = rep(0,TT+1)
mm[1]  = mm0
CC2[1] = CC20
for(i in 1:TT){
  mm[i+1]  = (y[i]/sigma2 + mm[i]/(tau2 + CC2[i]))/(1/sigma2 +
                                            1/(tau2 + CC2[i]))
  CC2[i+1] = 1/(1/sigma2 + 1/(tau2 + CC2[i]))
}


## Expected values
pdf(file.path("plots","Gaussian_means.pdf"))
plot(theta.t, type="n", xlab="t", ylab=expression(
                   paste("E(",theta[t],"|",y[1:t],")")), lwd=2)
lines(mm[-1], col="red", lwd=2, lty=2)
lines(mm.SIR[-1], col="blue", lwd=2, lty=3)
legend(120,3,c("Exact","SIR"), col=c("red","blue"), lty=c(2,3),
       lwd=2, bty="n")
dev.off()


## Confidence intervals
pdf(file.path("plots","Gaussian_confint.pdf"))
plot(theta.t, type="n", xlab="t", ylab=expression(
      paste("95% credible interval for ",theta[t])), lwd=2)
lines(mm[-1]-2*sqrt(CC2[-1]), col="red", lwd=1, lty=2)
lines(mm[-1]+2*sqrt(CC2[-1]), col="red", lwd=1, lty=2)
lines(mm.SIR[-1]-2*sqrt(CC2.SIR[-1]), col="blue", lwd=1, lty=3)
lines(mm.SIR[-1]+2*sqrt(CC2.SIR[-1]), col="blue", lwd=1, lty=3)
legend(120,3,c("Exact","SIR"), col=c("red","blue"), lty=c(2,3),
       lwd=2, bty="n")
dev.off()


## ESS as a function of time
pdf(file.path("plots","Gaussian_ESSintime.pdf"))
plot(ESS[-1], type="l", ylab="ESS", xlab="t",ylim=c(0,BB))
dev.off()


## ESS vs logarithm of the variance
pdf(file.path("plots","Gaussian_ESSvslogVar.pdf"))
plot(log(apply(weights.out,2,var)[-1]), ESS[-1], ylab="ESS",
     xlab=expression(paste("log(Var(",omega[t]^{(b)},"))")),
     ylim=c(0,BB))
dev.off()












