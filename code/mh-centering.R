## note can use this code to regenerate Abel's centering.pdf files


## generate simulated dataset
set.seed(1)
I = 10
J = 20
sigma2 = 1
tau2   = 2  

eta   = 0.8
theta = rnorm(I, eta, sqrt(tau2))
y     = matrix(NA, nrow=I, ncol=J)
for(i in 1:I){
  y[i,] = rnorm(J, theta[i], sqrt(sigma2))
}

nIts <- 10000

code <- nimbleCode({
    tau ~ dhalfflat()
    sigma ~ dhalfflat()
    for(i in 1:I) {
       for(j in 1:J) {
         y[i,j] ~ dnorm(eta + theta[i], sd = sigma)
         theta[i] ~ dnorm(0, sd = tau)
    }}
    eta ~ dflat()
})

lmeModel <- nimbleModel(code, data = list(y = y),
                 constants = list(I = I, J = J),
                 inits = list(tau = 1, sigma = 1, theta = rnorm(I), eta = 0))
clmeModel <- compileNimble(lmeModel)

conf <- configureMCMC(lmeModel, monitors = c('tau', 'sigma', 'eta', 'theta'))
conf$getSamplers()
## note: uniform on std deviation is equivalent to degenerate IG(-0.5, 0), hence conjugacy
mcmc <- buildMCMC(conf)
cmcmc <- compileNimble(mcmc, project = lmeModel)
set.seed(1)
cmcmc$run(nIts)
smpUncentered <- as.matrix(cmcmc$mvSamples)

pdf(file.path('plots','centering1.pdf'), width=4, height=4)
par(mfrow = c(2,2), mai=c(0.6,.5,.1,.1),mgp=c(1.8,.7,0))
ts.plot(smpUncentered[ , 'eta'], ylab = expression(eta))
ts.plot(smpUncentered[ , 'theta[1]'], ylab = expression(theta[1]))
acf(smpUncentered[ , 'eta'])
acf(smpUncentered[ , 'theta[1]'])
dev.off()

apply(smpUncentered, 2, coda::effectiveSize)

code <- nimbleCode({
    tau ~ dhalfflat()
    sigma ~ dhalfflat()
    for(i in 1:I) {
       for(j in 1:J) {
         y[i,j] ~ dnorm(theta[i], sd = sigma)
         theta[i] ~ dnorm(eta, sd = tau)
    }}
    eta ~ dflat()
})

lmeModel <- nimbleModel(code, data = list(y = y),
                 constants = list(I = I, J = J),
                 inits = list(tau = 1, sigma = 1, theta = rnorm(I), eta = 0))
clmeModel <- compileNimble(lmeModel)

conf <- configureMCMC(lmeModel, monitors = c('tau', 'sigma', 'eta', 'theta'))
conf$getSamplers()
## note: uniform on std deviation is equivalent to degenerate IG(-0.5, 0), hence conjugacy
mcmc <- buildMCMC(conf)
cmcmc <- compileNimble(mcmc, project = lmeModel)
set.seed(1)
cmcmc$run(nIts)
smpCentered <- as.matrix(cmcmc$mvSamples)

pdf(file.path('plots','centering2.pdf'), width=4, height=4)
par(mfrow = c(2,2), mai=c(0.6,.5,.1,.1),mgp=c(1.8,.7,0))
ts.plot(smpCentered[ , 'eta'], ylab = expression(eta))
ts.plot(smpCentered[ , 'theta[1]'] - smpCentered[ , 'eta'], ylab = expression(theta[1]))
acf(smpCentered[ , 'eta'])
acf(smpCentered[ , 'theta[1]'] - smpCentered[ , 'eta'])
dev.off()

apply(smpCentered, 2, coda::effectiveSize)


