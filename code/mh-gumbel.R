## @knitr gumbel-setup
n <- 100
x <- -log(-log(runif(n)))  # x is Gumbel-distributed with true theta = 0
### hist(x)

### Prior Parameters
xi <- 0
kappa2 <- 1

### Number of iterations
r <- 2000

### Initial value
thetaInit <- 0  # cheating a bit -- we start at the true value, so throwing away initial chain values won't be needed  

## @knitr gumbel-nimble-dist, warning=FALSE, message=FALSE
### Density of the Gumbel distribution, written as a nimbleFunction to allow use in BUGS code and compilation
dgumbel <- nimbleFunction(
  run = function(x = double(0), theta = double(0), log = integer(0, default = 0)) {
      returnType(double(0))
      logProb <- -(x-theta) - exp( -(x-theta))
      if(log){
          return(logProb)
      } else{
          return(exp(logProb))
      }
  })


## @knitr gumbel-bugs
gumbelModelCode <- nimbleCode({
    theta ~ dnorm(xi, sd = kappa)
    for(i in 1:n)
          x[i] ~ dgumbel(theta)
})   

gumbelModel <- nimbleModel(gumbelModelCode,
          constants = list(n = n, xi = xi, kappa = sqrt(kappa2)),
          data = list(x = x), inits = list(theta = thetaInit))

## @knitr gumbel-mcmc
### variance of proposal
tau2 <- 5
### empty MCMC configuration as default would use adaptive algorithm 
conf <- configureMCMC(gumbelModel, nodes = NULL)
### add basic Metropolis sampler for the parameter
conf$addSampler('theta', 'RW', control = list(adaptive = FALSE,
                                              scale = sqrt(tau2)))
### create MCMC algorithm for the model, given the configuration
mcmc <- buildMCMC(conf)
### compiled version of the model
cGumbelModel <- compileNimble(gumbelModel)
### compiled version of MCMC (model needs to be compiled first or at the same time)
cmcmc <- compileNimble(mcmc)
set.seed(1)
### run the MCMC
cmcmc$run(r)

## @knitr gumbel-results
### basic diagnostics
smp1 <- as.matrix(cmcmc$mvSamples)
## source('utils.R')
accRate1 <- calcAccRate(smp1[ , 1])
ess1 <- coda::effectiveSize(smp1)

## @knitr gumbel-other-prop-variances
tau2 <- 0.001
conf$removeSamplers()
conf$addSampler('theta', 'RW', control = list(adaptive = FALSE,
                                              scale = sqrt(tau2)))
mcmc <- buildMCMC(conf)
cmcmc <- compileNimble(mcmc, resetFunctions = TRUE)
cmcmc$run(r)

smp2 <- as.matrix(cmcmc$mvSamples)
accRate2 <- sum(diff(smp2) != 0) / r
ess2 <- coda::effectiveSize(smp2)

tau2 <- 0.07
conf$removeSamplers()
conf$addSampler('theta', 'RW', control = list(adaptive = FALSE, scale = sqrt(tau2)))
mcmc <- buildMCMC(conf)
cmcmc <- compileNimble(mcmc, resetFunctions = TRUE)
cmcmc$run(r)

smp3 <- as.matrix(cmcmc$mvSamples)
accRate3 <- sum(diff(smp3) != 0) / r
ess3 <- coda::effectiveSize(smp3)

## @knitr gumbel-plots
par(mfrow = c(1,3), cex.main = .7)
ts.plot(smp1[ , 1], ylab = expression(tau^2 == 5),
        main = paste0("Accept rate: ",
               round(accRate1, 2), "; ESS: ", round(ess1)))
ts.plot(smp2[ , 1], ylab = expression(tau^2 == 0.001),
        main = paste0("Accept rate: ",
               round(accRate2, 2), "; ESS: ", round(ess2)))
ts.plot(smp3[ , 1], ylab = expression(tau^2 == 0.07),
        main = paste0("Accept rate: ",
               round(accRate3, 2), "; ESS: ", round(ess3)))


## @knitr gumbel-density
dgumbel <- function(x, theta, log = FALSE){
  z <- -(x-theta) - exp( -(x-theta))
  if(log) return(z) else return(exp(z))
}
## @knitr gumbel-noNimble
smpR <- rep(0, r)
acceptRate <- 0
thetaCurr <- thetaInit
for(s in 1:r){
    thetaProp <- rnorm(1, thetaCurr, sqrt(tau2))
    rho <- sum(dgumbel(x, thetaProp, log = TRUE)) -
        sum(dgumbel(x, thetaCurr, log = TRUE))
    ## likelihood and prior evaluated on log scale to avoid underflow
    rho <- rho + dnorm(thetaProp, xi, sqrt(kappa2), log = TRUE) -
        dnorm(thetaCurr, xi, sqrt(kappa2), log = TRUE)
    u <- runif(1, 0, 1)
    if(log(u) < rho){  # compare on log scale and no need for min()
      thetaCurr <- thetaProp
      acceptRate <- acceptRate + 1
    }
    smpR[s] <- thetaCurr
}



## @knitr metropolis-nimble

sampler_RW <- nimbleFunction(
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
        ## control list extraction
        logScale      <- control$log
        scale         <- control$scale
        ## node list generation
        calcNodes  <- model$getDependencies(target)
    },
    run = function() {
        currentValue <- model[[target]]
        propLogScale <- 0
        if(logScale) {
            propLogScale <- rnorm(1, mean = 0, sd = scale)
            propValue <- currentValue * exp(propLogScale)
        } else propValue <- rnorm(1, mean = currentValue,  sd = scale)
        model[[target]] <<- propValue
        logMHR <- calculateDiff(model, calcNodes) + propLogScale
        jump <- decide(logMHR)
        if(jump) nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
        else     nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodes, logProb = TRUE)
    })



