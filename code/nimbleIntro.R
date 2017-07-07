## @knitr regression-example

regrCode <- nimbleCode({
    ## priors
    intcpt ~ dnorm(0, sd = 1000)
    slope ~ dnorm(0, sd = 1000)
    ## Gelman (2006) recommended prior:
    sigma ~ dunif(0, 100)  
    
    for(i in 1:4) {
        ## linear predictor (deterministic node)
        pred_y[i] <- intcpt + slope * x[i]
        ## likelihood 
        y[i] ~ dnorm(pred_y[i], sd = sigma)
    }
})

## @knitr litters-code
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
   }          })

## @knitr litters-model
consts <- list(G = 2, N = 16,
     n = matrix(c(13,12,12,11,9,10,9,9,8,11,8,10,13,10,12,9,
           10,9,10,5,9,9,13,7,5,10,7,6,10,10,10,7), nrow = 2))
data = list(r = matrix(c(13,12,12,11,9,10,9,9,8,10,8,9, 12,9,11,8,9,8,9,
                         4,8,7,11,4,4,5,5,3,7,3,7,0), nrow = 2))
inits <- list(a = c(2, 2), b = c(2, 2))
## create a model object that can be used and modified
model <- nimbleModel(glmmCode, constants = consts, data = data, inits = inits)

## @knitr litters-operating

model$p
model$simulate('p')
model$p
model$p <- matrix(0.5, consts$G, consts$N)
## IMPORTANT - don't assume you know what child nodes do or
## don't need to be updated once you've changed values in the model:
model$calculate(model$getDependencies('p'))
model$calculate('r')

## plot(model$graph)
model$getDependencies('a[1]')
model$getDependencies('a[1]', determOnly = TRUE)
model$getDependencies('a[1]', dataOnly = TRUE)
model$getDependencies('a[1]', dataOnly = TRUE, downstream = TRUE)

