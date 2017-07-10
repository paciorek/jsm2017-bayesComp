library(nimble)
### Hyperparameters
sigma2 = 1^2
tau2   = 2^2
mm0    = 0
CC20   = 4

### Nimble Gausian SSM Model Code
gaussianSSMCode <- nimbleCode({
  theta.t[1] ~ dnorm(mm0, sd = sqrt(CC20))
  for(i in 1:TT){
    theta.t[i+1] ~ dnorm(theta.t[i], sd = sqrt(tau2))
    y[i] ~ dnorm(theta.t[i+1], sd = sqrt(sigma2))
  }
})


### Generating data
TT     = 150
theta.t = rep(0, TT+1)
y       = rep(0, TT)
theta.t[1] = 2
for(i in 1:TT){
  theta.t[i+1] = rnorm(1, theta.t[i], sqrt(tau2))
  y[i] = rnorm(1, theta.t[i+1], sqrt(sigma2))
}

### Creating Nimble Gausian SSM 
gaussianSSM <- nimbleModel(code = gaussianSSMCode, data = list(y = y), 
                           constants = list(sigma2 = sigma2, tau2 = tau2, mm0 = mm0,
                                            CC20 = CC20, TT = 150),
                           dimensions = list(theta.t = 151, y = 150))
CgaussianSSM <- compileNimble(gaussianSSM)

### Creating SIR Filter.  "Thresh = 1" option ensures that resampling will happen at each time point.
gaussianSSM_SIRfilter <- buildBootstrapFilter(gaussianSSM, nodes = 'theta.t', control = list(thresh = 1, saveAll = TRUE))
CgaussianSSM_SIRfilter <- compileNimble(gaussianSSM_SIRfilter, project = gaussianSSM)
BB = 10000
CgaussianSSM_SIRfilter$run(BB)
weights.out <- exp(as.matrix(CgaussianSSM_SIRfilter$mvWSamples)[, 152:302])/apply(exp(as.matrix(CgaussianSSM_SIRfilter$mvWSamples)[, 152:302]), 2, sum)
particleSamples <- as.matrix(CgaussianSSM_SIRfilter$mvWSamples)[, 1:151]
mm.SIR <- apply(weights.out*particleSamples, 2, sum)
CC2.SIR <- apply(weights.out*particleSamples^2, 2, sum) - mm.SIR^2
ESS.SIR <- CgaussianSSM_SIRfilter$returnESS()

### Creating APF Filter.  
gaussianSSM_APfilter <- buildAuxiliaryFilter(gaussianSSM, nodes = 'theta.t', control = list(lookahead = 'mean', saveAll = TRUE))
CgaussianSSM_APfilter <- compileNimble(gaussianSSM_APfilter, project = gaussianSSM, resetFunctions = TRUE)
BB = 10000
CgaussianSSM_APfilter$run(BB)
weights.out <- exp(as.matrix(CgaussianSSM_APfilter$mvWSamples)[, 152:302])/apply(exp(as.matrix(CgaussianSSM_APfilter$mvWSamples)[, 152:302]), 2, sum)
particleSamples <- as.matrix(CgaussianSSM_APfilter$mvWSamples)[, 1:151]
mm.APF <- apply(weights.out*particleSamples, 2, sum)
CC2.APF <- apply(weights.out*particleSamples^2, 2, sum) - mm.APF^2
ESS.APF <- CgaussianSSM_APfilter$returnESS()

### Exact solution
mm  = rep(0,TT+1)
CC2 = rep(0,TT+1)
mm[1]  = mm0
CC2[1] = CC20
for(i in 1:TT){
  mm[i+1]  = (y[i]/sigma2 + mm[i]/(tau2 + CC2[i]))/(1/sigma2 + 1/(tau2 + CC2[i]))
  CC2[i+1] = 1/(1/sigma2 + 1/(tau2 + CC2[i]))
}

## Expected values
quartz()
plot(theta.t, type="n", xlab="t", ylab=expression(paste("E(",theta[t],"|",y[1:t],")")), lwd=2)
lines(mm[-1], col="red", lwd=2, lty=2)
lines(mm.SIR[-1], col="blue", lwd=2, lty=3)
lines(mm.APF[-1], col="green", lwd=2, lty=4)
legend(115,4,c("Exact","SIR","APF"), col=c("red","blue","green"), lty=c(2,3,4), lwd=2, bty="n")
dev.print(dev=pdf,file="Gaussian_means_APF_NIMBLE.pdf")


## Confidence intervals
quartz()
plot(theta.t, type="n", xlab="t", ylab=expression(paste("95% credible interval for ",theta[t])), lwd=2)
lines(mm[-1]-2*sqrt(CC2[-1]), col="red", lwd=1, lty=2)
lines(mm[-1]+2*sqrt(CC2[-1]), col="red", lwd=1, lty=2)
lines(mm.SIR[-1]-2*sqrt(CC2.SIR[-1]), col="blue", lwd=1, lty=3)
lines(mm.SIR[-1]+2*sqrt(CC2.SIR[-1]), col="blue", lwd=1, lty=3)
lines(mm.APF[-1]-2*sqrt(CC2.APF[-1]), col="green", lwd=1, lty=3)
lines(mm.APF[-1]+2*sqrt(CC2.APF[-1]), col="green", lwd=1, lty=3)
legend(115,4,c("Exact","SIR","APF"), col=c("red","blue","green"), lty=c(2,3,4), lwd=2, bty="n")
dev.print(dev=pdf,file="Gaussian_confint_APF_NIMBLE.pdf")


## ESS as a function of time
quartz()
plot(ESS.SIR, type="l", ylab="ESS", xlab="t", ylim=c(0,BB), lty=3, col="blue")
lines(ESS.APF,lty=4, col="green")
legend(120,9000,c("SIR","APF"), col=c("blue","green"), lty=c(3,4), lwd=2, bty="n")
dev.print(dev=pdf,file="Gaussian_ESSintime_APF_NIMBLE.pdf")

