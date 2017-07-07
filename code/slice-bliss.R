w <- c(1.6907, 1.7242, 1.7552, 1.7842, 1.8113, 1.8369, 1.8610, 1.8839)
y <- c(6, 13, 18, 28, 52, 53, 61, 60)
n <- c(59, 60, 62, 56, 63, 59, 62, 60)

hypers <- list(a0 = .25, b0 = .25, c0 = 2, d0 = 10, e0 = 2, f0 = .001)

dloggamma <- nimbleFunction(
    run = function(x = double(0), a = double(0), b = double(0),
        log = integer(0, default = 0)) {
        returnType(double(0))
        logProb <- a*x-b*exp(x)
        if(log) return(logProb)
        else return(exp(logProb)) 
        })

dlogsqrtIG <- nimbleFunction(
    run = function(x = double(0), a = double(0), b = double(0),
        log = integer(0, default = 0)) {
        returnType(double(0))
        logProb <- -2*a*x - b*exp(-2*x)
        if(log) return(logProb)
        else return(exp(logProb)) 
        })

blissCode <- nimbleCode({
    for(i in 1:G) {
        x[i] <- (w[i] - mu) / sigma
        p[i] <- (expit(x[i]))^m1
        y[i] ~ dbin(p[i], n[i])
    }
    m1 <- exp(logm1)
    logm1 ~ dloggamma(a0, b0)
    mu ~ dnorm(c0, sd = d0)
    logsigma ~ dlogsqrtIG(e0, f0)
    sigma <- exp(logsigma)
})

inits <- list(mu = 2, logsigma = -8, logm1 = -4)
# inits <- list(mu = -8, logsigma = 5, logm1 = 5)  ## bad starting values
blissModel <- nimbleModel(blissCode, data = list(w = w, y = y),
                 constants = c(list(n = n, G = length(n)), hypers),
                 inits = inits)
cBlissModel <- compileNimble(blissModel)

thin <- 1
nIts <- 5000
burnin <- 1000

## default (adaptive) samplers
conf <- configureMCMC(blissModel, thin = thin)
mcmc <- buildMCMC(conf)
cmcmc <- compileNimble(mcmc, project = blissModel)

set.seed(1)
cmcmc$run(nIts)

smp_default <- as.matrix(cmcmc$mvSamples)


conf <- configureMCMC(blissModel, nodes = NULL, thin = thin)
for(param in c('logm1', 'mu', 'logsigma'))
    conf$addSampler(param, 'slice')
conf$getSamplers()

mcmc <- buildMCMC(conf)
cmcmc <- compileNimble(mcmc, project = blissModel, resetFunctions = TRUE)
cBlissModel$setInits(inits)

set.seed(1)
cmcmc$run(nIts)

smp_slice <- as.matrix(cmcmc$mvSamples)

pdf(file.path('plots','slice-bliss.pdf'), width = 6, height = 3.5)
par(mfrow = c(2, 3), mai=c(0.4,.4,.3,.1),mgp=c(1.8,.7,0))
for(i in 1:3) {
    ts.plot(smp_default[ , i], ylab = colnames(smp_default)[i])
    if(i == 2) title('default univariate Metropolis')
}
for(i in 1:3) {
    ts.plot(smp_slice[ , i], ylab = colnames(smp_default)[i])
    if(i == 2) title('slice samplers')
}
 dev.off()          

smp_default <- smp_default[-seq_len(burnin), ]
cor(smp_default)
apply(smp_default, 2, coda::effectiveSize)

smp_slice <- smp_slice[-seq_len(burnin), ]
cor(smp_slice)
apply(smp_slice, 2, coda::effectiveSize)


