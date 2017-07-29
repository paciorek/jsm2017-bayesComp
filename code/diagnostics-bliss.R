## @knitr bliss-mcmc-twoStarts
nIts <- 3000

inits1 <- list(mu = 2, logsigma = -8, logm1 = -4)
inits2 <- list(mu = 1.7, logsigma = -1, logm1 = 0)
## inits3 <- list(mu = -8, logsigma = 5, logm1 = 5)  ## bad starting values - no convergence after 40k iterations if using these

c2 <- (2.4/sqrt(length(params)))^2 
propCov <- c2 * cov(smpSave)

## chain 1
conf$getSamplers()
conf$removeSamplers()
conf$addSampler(params, type = 'RW_block',
                control = list(adaptive = FALSE, scale = 1,
                               propCov = propCov))

## rebuild and recompile MCMC as we have changed samplers
mcmc <- buildMCMC(conf)
cmcmc <- compileNimble(mcmc, project = blissModel, resetFunctions = TRUE)

cBlissModel$setInits(inits1)
set.seed(1)
cmcmc$run(nIts)
smp1 <- as.matrix(cmcmc$mvSamples)

cBlissModel$setInits(inits2)
set.seed(1)
cmcmc$run(nIts)

smp2 <- as.matrix(cmcmc$mvSamples)

## note that NIMBLE can also run multiple MCMCs using runMCMC

## @knitr bliss-gelman-rubin

library(coda, warn.conflicts = FALSE)

chains <- as.mcmc.list(
    list(as.mcmc(smp1), as.mcmc(smp2)))
gr = gelman.diag(chains, autoburnin = FALSE)
print(gr)

## @knitr bliss-gelman-rubin-plot
pdf(file.path('plots', 'bliss-gelman-rubin.pdf'), width=5, height=4)
par(mfrow=c(3, 3), mai=c(0.4,.4,.3,.1),mgp=c(1.8,.7,0))
gelman.plot(chains, autoburnin = FALSE, auto.layout = FALSE,
            col = 1:3)
nm <- colnames(smp1)
for(i in 1:3)
    ts.plot(smp1[, i], main = nm[i], xlab = '', ylab = '')
for(i in 1:3)
    ts.plot(smp2[, i], main = nm[i], xlab = '', ylab = '')
dev.off()

## chain not coverged during initial few hundred samples, so R value is still large
## let's see what happens if we exclude burnin (by default half the chain)

gr = gelman.diag(chains, autoburnin = TRUE)
print(gr)
## force plot to only assess last half of chains:
chains <- as.mcmc.list(
    list(as.mcmc(smp1[1501:3000, ]), as.mcmc(smp2[1501:3000, ])))
gelman.plot(chains, autoburnin = FALSE, auto.layout = FALSE,
            col = 1:3)
