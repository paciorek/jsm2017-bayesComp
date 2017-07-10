library(nimble)

### Hyperparameters
eta    = log(0.15/sqrt(252))
rho    = 0.95
tau2   = 0.25


### Nimble Stochastic Volatility Model Code
svCode <- nimbleCode({
  theta.t[1] ~ dnorm(eta, sd = sqrt(tau2/(1-rho^2)))
  for(i in 1:TT){
    theta.t[i+1] ~ dnorm(eta + rho*(theta.t[i] - eta), sd = sqrt(tau2))
    y[i] ~ dnorm(0, sd = sqrt(exp(2*theta.t[i+1])))
  }
})

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

### Creating Nimble SV model 
svModel <- nimbleModel(code = svCode, data = list(y = y), 
                           constants = list(eta = eta, rho = rho, tau2 = tau2, TT = 150),
                           dimensions = list(theta.t = 151, y = 150))
CsvModel <- compileNimble(svModel)

### Creating SIR Filter.  "Thresh = 1" option ensures that resampling will happen at each time point.
svModel_SIRfilter <- buildBootstrapFilter(svModel, nodes = 'theta.t', control = list(thresh = 1, saveAll = TRUE))
CsvModel_SIRfilter <- compileNimble(svModel_SIRfilter, project = svModel)
BB = 10000
CsvModel_SIRfilter$run(BB)
weights.out <- exp(as.matrix(CsvModel_SIRfilter$mvWSamples)[, 152:302])/apply(exp(as.matrix(CsvModel_SIRfilter$mvWSamples)[, 152:302]), 2, sum)
particleSamples <- as.matrix(CsvModel_SIRfilter$mvWSamples)[, 1:151]
mm.SIR <- apply(weights.out*particleSamples, 2, sum)
CC2.SIR <- apply(weights.out*particleSamples^2, 2, sum) - mm.SIR^2
ESS.SIR <- CsvModel_SIRfilter$returnESS()

### Creating APF Filter.  
svModel_APfilter <- buildAuxiliaryFilter(svModel, nodes = 'theta.t', control = list(lookahead = 'mean', saveAll = TRUE))
CsvModel_APfilter <- compileNimble(svModel_APfilter, project = svModel, resetFunctions = TRUE)
BB = 10000
CsvModel_APfilter$run(BB)
weights.out <- exp(as.matrix(CsvModel_APfilter$mvWSamples)[, 152:302])/apply(exp(as.matrix(CsvModel_APfilter$mvWSamples)[, 152:302]), 2, sum)
particleSamples <- as.matrix(CsvModel_APfilter$mvWSamples)[, 1:151]
mm.APF <- apply(weights.out*particleSamples, 2, sum)
CC2.APF <- apply(weights.out*particleSamples^2, 2, sum) - mm.APF^2
ESS.APF <- CsvModel_APfilter$returnESS()

## Expected values
quartz()
plot(theta.t, type="l", xlab="t", ylab=expression(paste("E(",theta[t],"|",y[1:t],")")), lwd=2)
lines(mm.SIR[-1], col="blue", lwd=2, lty=3)
lines(mm.APF[-1], col="green", lwd=2, lty=4)
legend(120,4,c("True","SIR","APF"), col=c("black","blue","green"), lty=c(1,3,4), lwd=2, bty="n")
#dev.print(dev=pdf,file="Gaussian_means_APF_NIMBLE.pdf")


## ESS as a function of time
quartz()
plot(ESS.SIR, type="l", ylab="ESS", xlab="t", ylim=c(0,BB), lty=3, col="blue")
lines(ESS.APF,lty=4, col="green")
legend(120,9000,c("SIR","APF"), col=c("blue","green"), lty=c(3,4), lwd=2, bty="n")
#dev.print(dev=pdf,file="Gaussian_ESSintime_APF_NIMBLE.pdf")

mean(ESS.SIR)
mean(ESS.APF)












