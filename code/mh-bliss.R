## @knitr bliss-setup
## data 
w <- c(1.6907, 1.7242, 1.7552, 1.7842, 1.8113, 1.8369, 1.8610, 1.8839)
y <- c(6, 13, 18, 28, 52, 53, 61, 60)
n <- c(59, 60, 62, 56, 63, 59, 62, 60)

## plot(w, y/n)

hypers <- list(a0 = .25, b0 = .25, c0 = 2, d0 = 10, e0 = 2, f0 = .001)

## user-defined distributions for the transformed priors on
## log sigma and log m1

## implied prior on log sigma
dlogsqrtIG <- nimbleFunction(
    run = function(x = double(0), a = double(0), b = double(0),
        log = integer(0, default = 0)) {
        returnType(double(0))
        logProb <- -2*a*x - b*exp(-2*x)
        if(log) return(logProb)
        else return(exp(logProb)) 
        })

## implied prior on log m1
dloggamma <- nimbleFunction(
    run = function(x = double(0), a = double(0), b = double(0),
        log = integer(0, default = 0)) {
        returnType(double(0))
        logProb <- a*x-b*exp(x)
        if(log) return(logProb)
        else return(exp(logProb)) 
    })

## @knitr bliss-model

## setup BUGS code
blissCode <- nimbleCode({
    mu ~ dnorm(c0, sd = d0)
    logsigma ~ dlogsqrtIG(e0, f0)
    logm1 ~ dloggamma(a0, b0)
    sigma <- exp(logsigma)
    m1 <- exp(logm1)
    for(i in 1:G) {
        x[i] <- (w[i] - mu) / sigma
        p[i] <- (expit(x[i]))^m1
        y[i] ~ dbin(p[i], n[i])
    }
})

## create model structures
inits <- list(mu = 2, logsigma = -8, logm1 = -4) 
# inits <- list(mu = -8, logsigma = 5, logm1 = 5)  ## bad starting values
blissModel <- nimbleModel(blissCode, data = list(w = w, y = y),
                 constants = c(list(n = n, G = length(n)), hypers),
                 inits = inits)
cBlissModel <- compileNimble(blissModel)

## @knitr bliss-mcmc-trialrun
## run MCMC with diagonal covariance to estimate posterior covariance
params <- blissModel$getNodeNames(stochOnly = TRUE, includeData = FALSE)
conf <- configureMCMC(blissModel)
conf$removeSamplers()  ## remove NIMBLE's default samplers
conf$addSampler(params, type = 'RW_block',
                control = list(adaptive = FALSE, scale = 1,
                               propCov = diag(rep(0.01, 3))))
conf$getSamplers()

mcmc <- buildMCMC(conf)
cmcmc <- compileNimble(mcmc, project = blissModel)

nIts <- 10000
set.seed(1)
cmcmc$run(nIts)

## examine posterior samples
smp <- as.matrix(cmcmc$mvSamples)[ , params]
par(mfrow = c(1, 3))
for(i in 1:3)
      ts.plot(smp[ , i])

burnin <- 1000
smpSave <- smp[-seq_len(burnin), ]
cor(smpSave)
apply(smpSave, 2, coda::effectiveSize)

## @knitr bliss-mcmc-theoretical
## then use blocked Metropolis with theoretical proposal covariance

c2 <- (2.4/sqrt(length(params)))^2 
theoreticalPropCov <- c2 * cov(smpSave)

conf$getSamplers()
conf$removeSamplers()
conf$addSampler(params, type = 'RW_block',
                control = list(adaptive = FALSE, scale = 1,
                               propCov = theoreticalPropCov))

mcmc <- buildMCMC(conf)
cmcmc <- compileNimble(mcmc, project = blissModel, resetFunctions = TRUE)

nIts <- 9000
set.seed(1)
cmcmc$run(nIts)

smp <- as.matrix(cmcmc$mvSamples)
apply(smp, 2, calcAccRate)

par(mfrow = c(1, 3))
for(i in 1:3)
      ts.plot(smp[ , i])

apply(smp, 2, coda::effectiveSize)

## now illustrate adaptive Metropolis by comparing adaptive univariate,
## adaptive blocked, and theoretical blocked samplers

nIts <- 10000

## @knitr bliss-mcmc-blocked-theory
## blocked Metropolis with theoretical proposal covariance

conf <- configureMCMC(blissModel, nodes = NULL)
conf$addSampler(params, type = 'RW_block',
                control = list(adaptive = FALSE, scale = 1,
                               propCov = theoreticalPropCov))

mcmc <- buildMCMC(conf)
cmcmc <- compileNimble(mcmc, project = blissModel, resetFunctions = TRUE)

set.seed(1)
cBlissModel$setInits(inits)
cmcmc$run(nIts)

smp1 <- as.matrix(cmcmc$mvSamples)[ , params]
apply(smp, 2, calcAccRate)

par(mfrow = c(1, 3))
for(i in 1:3)
      ts.plot(smp1[ , i])

burnin <- 1000
apply(smp1[-seq_len(burnin), ], 2, coda::effectiveSize)

## @knitr bliss-mcmc-blocked-adaptive
## blocked Metropolis with adaptive proposal covariance

## this will allow us to see the proposal covariance the adaptive
## scheme evolves to
nimbleOptions(buildInterfacesForCompiledNestedNimbleFunctions = TRUE)

## now we set up adaptive block Metropolis-Hastings
## sometimes adaptive block sampler performs poorly if initial
## proposal scales are poor so set them here, but
## this is not always necessary and should be largely unnecessary
## as of fall 2017 when we will put a new scheme in place

conf$removeSamplers()
conf$addSampler(params, type = 'RW_block',
                control = list(adaptive = TRUE, scale = 1,
                               propCov = diag(rep(0.001, 3))))

mcmc <- buildMCMC(conf)
cmcmc <- compileNimble(mcmc, project = blissModel, resetFunctions = TRUE)

set.seed(1)
cBlissModel$setInits(inits)
cmcmc$run(nIts)
## adapted proposal covariance
cov2cor(cmcmc$samplerFunctions[[1]]$propCov)

smp2 <- as.matrix(cmcmc$mvSamples)[ , params]
for(i in 1:3)
      ts.plot(smp2[ , i])

cor(smp2[5000:9000, ])
burnin <- 5000
apply(smp2[-seq_len(burnin), ], 2, coda::effectiveSize)

## @knitr bliss-mcmc-univariate-adaptive
## univariate Metropolis with adaptive proposal scale (sd)
conf <- configureMCMC(blissModel) ## univariate adaptive is default
conf$getSamplers()

mcmc <- buildMCMC(conf)
cmcmc <- compileNimble(mcmc, project = blissModel, resetFunctions = TRUE)

set.seed(1)
cBlissModel$setInits(inits)
cmcmc$run(nIts)

smp3 <- as.matrix(cmcmc$mvSamples)[ , params]
par(mfrow = c(1, 3))
for(i in 1:3)
      ts.plot(smp3[ , i])

burnin <- 1000
apply(smp3[-seq_len(burnin), ], 2, coda::effectiveSize)

## @knitr bliss-mcmc-comparison-plots

pdf(file.path('plots','adaptive-bliss.pdf'), width=5, height=3.5)
par(mfrow = c(3,3), mai=c(0.35,.4,.2,.1),mgp=c(1.8,.7,0), cex.main = 1)
for(i in 1:3) {
    ts.plot(smp1[ , i], ylab = params[i])
    if(i == 2) title('theoretical covariance')
}
for(i in 1:3) {
    ts.plot(smp2[ , i], ylab = params[i])
    if(i == 2) title('adaptive block')
}
for(i in 1:3) {
    ts.plot(smp3[ , i], ylab = params[i])
    if(i == 2) title('adaptive univariate')
}
dev.off()
