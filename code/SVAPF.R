library(nimble)

### Hyperparameters
eta    = log(0.15/sqrt(252))
rho    = 0.95
tau2   = 0.25


### Nimble Stochastic Volatility Model Code
sv.Code <- nimbleCode({
  theta.t[1] ~ dnorm(eta, sd = sqrt(tau2/(1-rho^2)))
  for(i in 1:TT){
    theta.t[i+1] ~ dnorm(eta + rho*(theta.t[i] - eta), sd = sqrt(tau2))
    y[i] ~ dnorm(0, sd = sqrt(exp(2*theta.t[i+1])))
  }
})

### Generating data
set.seed(1)
TT     = 150
theta.t = rep(0, TT+1)
y       = rep(0, TT)

theta.t[1] = eta
for(i in 1:TT){
  theta.t[i+1] = rnorm(1, eta + rho*(theta.t[i] - eta), sqrt(tau2))
  y[i] = rnorm(1, 0, exp(theta.t[i+1]))
}

plot(y, type="l")

### Creating Nimble SV model 
sv.Model <- nimbleModel(code = sv.Code, data = list(y = y), 
                           constants = list(eta = eta, rho = rho, 
                                            tau2 = tau2, TT = 150),
                           dimensions = list(theta.t = 151, y = 150))
Csv.Model <- compileNimble(sv.Model)

### Creating SIR Filter.  
### "thresh = 1" option ensures that resampling will happen at each time point.
sv.Model.SIRfilter <- buildBootstrapFilter(sv.Model, nodes = 'theta.t',
                                           control = list(thresh = 1, 
                                                          saveAll = TRUE))
Csv.Model.SIRfilter <- compileNimble(sv.Model.SIRfilter, project = sv.Model)
BB = 10000
Csv.Model.SIRfilter$run(BB)

### The first 151 columns of the output will be particle samples 
### for all 151 time points.
particle.samples <- as.matrix(Csv.Model.SIRfilter$mvWSamples)[, 1:151]

### Columns 152 through 302 of the output will be the unnormalized log weights
### for each particle at each time point.
log.weights.out <- as.matrix(Csv.Model.SIRfilter$mvWSamples)[, 152:302]
### Below, we un-log and normalize the weights
weights.out <-  exp(log.weights.out)/apply(exp(log.weights.out), 2, sum)

mm.SIR <- apply(weights.out*particle.samples, 2, sum)
CC2.SIR <- apply(weights.out*particle.samples^2, 2, sum) - mm.SIR^2
ESS.SIR <- Csv.Model.SIRfilter$returnESS()

### ESS averaged over all time points. Time point 1 corresponds to the 
### prior distribution placed on theta in the model, so we drop it
print(mean(ESS.SIR[-1]))

### Creating APF Filter.  
svModel.APfilter <- buildAuxiliaryFilter(sv.Model, nodes = 'theta.t',
                                         control = list(lookahead = 'mean',
                                                        saveAll = TRUE))
CsvModel.APfilter <- compileNimble(svModel.APfilter, project = sv.Model,
                                   resetFunctions = TRUE)
BB = 10000
CsvModel.APfilter$run(BB)
particle.samples <- as.matrix(CsvModel.APfilter$mvWSamples)[, 1:151]
log.weights.out <- as.matrix(CsvModel.APfilter$mvWSamples)[, 152:302]
weights.out <-  exp(log.weights.out)/apply(exp(log.weights.out), 2, sum)

mm.APF <- apply(weights.out*particle.samples, 2, sum)
CC2.APF <- apply(weights.out*particle.samples^2, 2, sum) - mm.APF^2
ESS.APF <- CsvModel.APfilter$returnESS()

print(mean(ESS.APF[-1]))

## Expected values
plot(theta.t, type="l", xlab="t", 
     ylab=expression(paste("E(",theta[t],"|",y[1:t],")")), lwd=2)
lines(mm.SIR[-1], col="blue", lwd=2, lty=3)
lines(mm.APF[-1], col="green", lwd=2, lty=4)
legend(120,-2,c("True","SIR","APF"), col=c("black","blue","green"),
       lty=c(1,3,4), lwd=2, bty="n")


## ESS as a function of time
plot(ESS.SIR, type="l", ylab="ESS", xlab="t", ylim=c(0,BB), lty=3, col="blue")
lines(ESS.APF,lty=4, col="green")
legend(120,10000,c("SIR","APF"), col=c("blue","green"), lty=c(3,4),
       lwd=2, bty="n")













