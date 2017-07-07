## @knitr litters-code2
glmmCode <- nimbleCode({
  for (i in 1:G) {
      ## priors for hyperparameters (from original example, not necessarily recommended)
     a[i] ~ dgamma(1, .001)
     b[i] ~ dgamma(1, .001)
     for (j in 1:N) {
	## random effects ('latent process')
        p[i,j] ~ dbeta(a[i], b[i]) 
     	## likelihood 
        r[i,j] ~ dbin(p[i,j], n[i,j])
     }
  }
})

## @knitr litters-model2
consts <- list(G = 2, N = 16,
     n = matrix(c(13,12,12,11,9,10,9,9,8,11,8,10,13,10,12,9,
           10,9,10,5,9,9,13,7,5,10,7,6,10,10,10,7), nrow = 2))
data = list(r = matrix(c(13,12,12,11,9,10,9,9,8,10,8,9, 12,9,11,8,9,8,9,
                         4,8,7,11,4,4,5,5,3,7,3,7,0), nrow = 2))
inits <- list(a = c(2, 2),
              b = c(2, 2))
## create a model object that can be used and modified
glmmModel <- nimbleModel(glmmCode, constants = consts, data = data, inits = inits)
cglmmModel <- compileNimble(glmmModel)

## @knitr litters-default-mcmc

thin <- 10
# thinning _only_ to reduce time of plotting 
conf <- configureMCMC(glmmModel, thin = thin)
conf$getSamplers(c('a','b','p[1,1]','p[2,1]'))

## @knitr litters-run-mcmc

conf$addMonitors(c('a', 'b', 'p'))
mcmc <- buildMCMC(conf)

cmcmc <- compileNimble(mcmc, project = glmmModel)

nIts <- 10000
set.seed(1)
cmcmc$run(nIts)

smp_basic <- as.matrix(cmcmc$mvSamples)

## @knitr litters-mcmc-results
pdf(file.path('plots','gibbs-litters.pdf'), width=6, height=4)
par(mfrow = c(2, 4), mai=c(0.3,.2,.4,.1),mgp=c(1.8,.7,0))
ts.plot(smp_basic[ , 'a[1]'], main = expression(a[1]), xlab = '', ylab = '')
ts.plot(smp_basic[ , 'b[1]'], main = expression(b[1]), xlab = '', ylab = '')
ts.plot(smp_basic[ , 'a[2]'], main = expression(a[2]), xlab = '', ylab = '', ylim = c(0,20))
ts.plot(smp_basic[ , 'b[2]'], main = expression(b[2]), xlab = '', ylab = '', ylim = c(0,6))
ts.plot(smp_basic[ , 'p[1, 3]'], main = "p[1,3]", xlab = '', ylab = '', ylim = c(0,1))
ts.plot(smp_basic[ , 'p[1, 15]'], main = "p[1,15]", xlab = '', ylab = '', ylim = c(0,1))
ts.plot(smp_basic[ , 'p[2, 3]'], main = "p[2,3]", xlab = '', ylab = '', ylim = c(0,1))
ts.plot(smp_basic[ , 'p[2, 15]'], main = "p[2,15]", xlab = '', ylab = '', ylim = c(0,1))
dev.off()
