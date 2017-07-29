## @knitr faithful-data

y <- faithful$waiting
n <- length(y)
## hist(y)

## @knitr mixture-auxiliary-code

## use of non-constant indexes requires NIMBLE version 0.6-6
if(packageVersion('nimble') < '0.6.6')
    stop("NIMBLE version 0.6-6 or greater required for non-constant indexes.")

nimbleOptions(allowDynamicIndexing = TRUE) ## it's a beta feature in 0.6-6

mixCode <- nimbleCode({
    for(i in 1:n) {
        y[i] ~ dnorm(theta[ksi[i]+1], var = sigma2[ksi[i]+1])
        ksi[i] ~ dbern(omega)
    }
    omega ~ dbeta(1, 1)
    for(j in 1:2) {
        theta[j] ~ dnorm(mu, tau2)
        sigma2[j] ~ dinvgamma(a, c)
    }
})

## @knitr fit-mixture-auxiliary

mixModel <- nimbleModel(mixCode, data = list(y = y),
                        constants = list(n = n, mu = 0, tau2 = .000001,
                                         a = .001, c = .001),
                     inits = list(omega = 0.5, theta = rep(mean(y), 2),
                                  sigma2 = rep(var(y), 2),
                                  ksi = sample(c(0,1), n, replace = T)))

cmixModel <- compileNimble(mixModel)

conf <- configureMCMC(mixModel, monitors = c('ksi', 'omega', 'theta', 'sigma2'))
mcmc <- buildMCMC(conf)
cmcmc <- compileNimble(mcmc, project = mixModel)

set.seed(1)
cmcmc$run(1000)

smp <- as.matrix(cmcmc$mvSamples)

burnin <- 500
smp <- smp[-seq_len(burnin), ]
postMeans <- colMeans(smp)
postMeans
mean(y > 65)

xgrid <- seq(min(y), max(y), len = 100)
predDens <- (1-postMeans['omega'])*dnorm(xgrid, postMeans['theta[1]'], postMeans['sigma[1]']) +
         postMeans['omega']*dnorm(xgrid, postMeans['theta[2]'], postMeans['sigma[2]'])

hist(y, freq = FALSE)
lines(xgrid, predDens)

pdf(file.path('plots','mixture-faithful.pdf'), width=6, height=2)
par(mfrow = c(1,4), mai=c(0.5,.5,.1,.1), mgp=c(1.8,.7,0),
    cex.main = 0.7)
hist(y, xlab = 'data')
ts.plot(smp[ , 'omega'], ylab = expression(omega),
        main = 'mixture weight')
ts.plot(1 + smp[ , 'ksi[201]'], ylab = expression(ksi[y==60]),
        main = 'mixture membership, y=60')
ts.plot(1 + smp[ , 'ksi[33]'], ylab = expression(ksi[y==66]),
        main = 'mixture membership, y=66')
dev.off()

## @knitr fit-mixture-marginalized

## user-defined two-component normal mixture density
dmix2N <- nimbleFunction(run = function(x = double(1),
                                       omega = double(0),
                                       theta = double(1), sigma = double(1),
                                       log = integer(0, default = 0)) {
    returnType(double(0))
    K <- 2
    omegaFull <- numeric(K)
    omegaFull[1] <- omega
    omegaFull[2] <- 1-omega
    n <- length(x)
    logdens <- 0
    for(i in 1:n) {
        tmpdens <- 0
        for(k in 1:K) {
            tmpdens <- tmpdens + omegaFull[k] * dnorm(x[i], theta[k], sigma[k], log = FALSE)
        }
        logdens <- logdens + log(tmpdens)
    }
        
    if(log) return(logdens) else return(exp(logdens))
})

mixCode <- nimbleCode({
    y[1:n] ~ dmix2N(omega, theta[1:2], sigma[1:2])
    omega ~ dunif(0, 1)
    for(k in 1:2) {
        theta[k] ~ dnorm(0, sd = 1000)
        sigma[k] ~ dhalfflat()
    }
})

mixModel <- nimbleModel(mixCode, data = list(y = y),
                     constants = list(n = n),
                     inits = list(omega = 0.5, theta = rep(mean(y), 2),
                                  sigma = rep(sd(y), 2)))

cmixMmodel <- compileNimble(mixModel)

## can build MCMC without configuration if ok with default samplers
mcmc <- buildMCMC(mixModel, monitors = c('omega', 'theta', 'sigma'))
cmcmc <- compileNimble(mcmc, project = mixModel)

set.seed(1)
cmcmc$run(5000)

smp <- as.matrix(cmcmc$mvSamples)

burnin <- 1000
smp <- smp[-seq_len(burnin), ]
postMeans <- colMeans(smp)
postMeans
mean(y > 65)

xgrid <- seq(min(y), max(y), len = 100)
predDens <- postMeans['omega']*dnorm(xgrid, postMeans['theta[1]'], postMeans['sigma[1]']) +
         (1-postMeans['omega'])*dnorm(xgrid, postMeans['theta[2]'], postMeans['sigma[2]'])

hist(y, freq = FALSE)
lines(xgrid, predDens)

