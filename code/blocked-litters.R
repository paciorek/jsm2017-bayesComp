## @knitr litters-blocked-config
thin <- 10

conf <- configureMCMC(glmmModel, nodes = NULL,
                      monitors = c('a','b','p'), thin = thin)
abNodes <- glmmModel$expandNodeNames(c('a', 'b'))
conf$addSampler(c('a[1]', 'b[1]'), 'crossLevel')
conf$addSampler(c('a[2]', 'b[2]'), 'crossLevel')

## @knitr litters-blocked-mcmc
nIts <- 10000

mcmc <- buildMCMC(conf)
cmcmc <- compileNimble(mcmc, project = glmmModel, resetFunctions = TRUE)

cglmmModel$setInits(list(a = c(2, 2), b=c(2, 2)))
set.seed(1)
cmcmc$run(nIts)

smp_blocked <- as.matrix(cmcmc$mvSamples)

pdf(file.path('plots', 'blocked-litters.pdf'), width=6, height=3)
par(mfrow = c(2,4), mai=c(0.3,.2,.4,.1),mgp=c(1.8,.7,0))
ts.plot(smp_blocked[ , 'a[1]'], main = expression(a[1]), xlab = '', ylab = '')
ts.plot(smp_blocked[ , 'b[1]'], main = expression(b[1]), xlab = '', ylab = '')
ts.plot(smp_blocked[ , 'a[2]'], main = expression(a[2]), xlab = '', ylab = '', ylim = c(0,20))
ts.plot(smp_blocked[ , 'b[2]'], main = expression(b[2]), xlab = '', ylab = '', ylim = c(0,6))
ts.plot(smp_blocked[ , 'p[1, 3]'], main = "p[1,3]", xlab = '', ylab = '', ylim = c(0,1))
ts.plot(smp_blocked[ , 'p[1, 15]'], main = "p[1,15]", xlab = '', ylab = '', ylim = c(0,1))
ts.plot(smp_blocked[ , 'p[2, 3]'], main = "p[2,3]", xlab = '', ylab = '', ylim = c(0,1))
ts.plot(smp_blocked[ , 'p[2, 15]'], main = "p[2,15]", xlab = '', ylab = '', ylim = c(0,1))
dev.off()
