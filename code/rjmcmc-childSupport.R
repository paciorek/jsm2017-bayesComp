source('rjmcmc-sampler.R')
source('rjmcmc-childSupport-generateData.R')

## sanity check against a GLM that ignores city effect (no city random effect)
personsGLM = persons
personsGLM$enforce <- cities$enforce[persons$city_ID]
personsGLM$benefit <- cities$benefit[persons$city_ID]
summary(glm(support ~ dad_age + enforce + benefit, data = personsGLM,
            family=binomial(link='probit')))


codeHierarchical <- nimbleCode({
    beta_benefit ~ dnorm(0, sd = 100)
    beta_enforce ~ dnorm(0, sd = 3) # model selection sensitive to vague priors
    beta_age ~ dnorm(0, sd = 100)
    beta0 ~ dnorm(0, sd = 100)
    sigma_alpha ~ dhalfflat()
    inc ~ dbern(0.5)  ## indicator variable for including the variable of interest
    betaInc <- beta_enforce * inc
    for(i in 1:n) {
        prob[i] <- pnorm(beta0 + beta_age * dad_age[i] + alpha[j[i]])
        support[i] ~ dbern(prob[i])
    }
    ## city-level model, including random effects
    for(k in 1:m) 
        alpha[k] ~ dnorm(betaInc * enforce[k] + beta_benefit * benefit[k],
                         sd = sigma_alpha)
})

## this parameterization, with the hierarchical regression terms (city-level regressors)
## moved to the data level, gives better mixing even without data augmentation 
code <- nimbleCode({
    beta_benefit ~ dnorm(0, sd = 100)
    beta_enforce ~ dnorm(0, sd = 3) # model selection sensitive to vague priors
    beta_age ~ dnorm(0, sd = 100)
    beta0 ~ dnorm(0, sd = 100)
    sigma_alpha ~ dhalfflat()
    inc ~ dbern(0.5)  ## indicator variable for including the variable of interest
    betaInc <- beta_enforce * inc
    for(i in 1:n) {
        ## city-level terms in the main likelihood
        prob[i] <- pnorm(beta0 + beta_age * dad_age[i] + betaInc * enforce[j[i]] +
            beta_benefit * benefit[j[i]] + alpha[j[i]])
        support[i] ~ dbern(prob[i])
    }
    ## city random effects
    for(k in 1:m) 
        alpha[k] ~ dnorm(0, sd = sigma_alpha)
})

model <- nimbleModel(code, constants = list(n = n, m = m, dad_age = persons$dad_age,
                                        enforce = cities$enforce,
                                        benefit  = cities$benefit,
                                        j = persons$city_ID),
                 data = list(support = persons$support),
                 inits = list(beta_benefit = 0, beta_enforce = 0, beta_age = 0,
                              beta0 = 0, sigma_alpha = 1, inc = 1))
cmodel <- compileNimble(model)

conf <- configureMCMC(model)
conf$getSamplers()

conf$removeSamplers('inc')
conf$addSampler(target = 'inc', type = RJindicatorSampler,
              control = list(scale = 1, coef = 'beta_enforce'))
conf$removeSamplers('beta_enforce')
conf$addSampler(target = 'beta_enforce', type = 'RW_sampler_nonzero_indicator',
                control = list(indicator = 'inc',
                           control = list(adaptive = TRUE, adaptInterval = 100,
                                  scale = .5, log = FALSE, reflective = FALSE)))

mcmc <- buildMCMC(conf)
cmcmc <- compileNimble(mcmc, project = model)

cmcmc$run(5000)

smp <- as.matrix(cmcmc$mvSamples)

beta_enforce <- smp[ , 'beta_enforce']

par(mfrow = c(2,3))
ts.plot(smp[ , 'beta_benefit'])
ts.plot(smp[ , 'beta_age'])
ts.plot(smp[ , 'sigma_alpha'])
ts.plot(beta_enforce)
## beta_enforce when it is in the model
ts.plot(beta_enforce[beta_enforce != 0])  
hist(beta_enforce[beta_enforce != 0])

burnin <- 1000
smp <- smp[(burnin+1):nrow(smp), ]
## inclusion probability for enforcement
mean(smp[ , 'inc'])

## now try Albert&Chib data augmentation trick

codeDataAug <- nimbleCode({
    beta_benefit ~ dnorm(0, sd = 100)
    beta_enforce ~ dnorm(0, sd = 3) # model selection sensitive to vague priors
    beta_age ~ dnorm(0, sd = 100)
    beta0 ~ dnorm(0, sd = 100)
    sigma_alpha ~ dhalfflat()
    inc ~ dbern(0.5)  ## indicator variable for including the variable of interest
    betaInc <- beta_enforce * inc
    for(i in 1:n) {
        ## Albert & Chib probit regression data augmentation trick
        z[i] ~ dnorm(beta0 + beta_age * dad_age[i] + alpha[j[i]], 1)
        support[i] ~ dinterval(z[i], 0)  # if support[i] = 0, z[i] < 0
    }
    for(k in 1:m) 
        alpha[k] ~ dnorm(betaInc * enforce[k] + beta_benefit * benefit[k],
                         sd = sigma_alpha)
})

model <- nimbleModel(codeDataAug, constants = list(n = n, m = m, dad_age = persons$dad_age,
                                        enforce = cities$enforce,
                                        benefit  = cities$benefit,
                                        j = persons$city_ID),
                 data = list(support = persons$support),
                 inits = list(beta_benefit = 0, beta_enforce = 0, beta_age = 0,
                              beta0 = 0, sigma_alpha = 1, inc = 1,
                              z = ifelse(persons$support == 0, -0.5, .5)))
cmodel <- compileNimble(model)

conf <- configureMCMC(model)
conf$getSamplers()[1:30]

conf$removeSamplers('inc')
conf$addSampler(target = 'inc', type = RJindicatorSampler,
              control = list(scale = 1, coef = 'beta_enforce'))
conf$removeSamplers('beta_enforce')
conf$addSampler(target = 'beta_enforce', type = 'RW_sampler_nonzero_indicator',
                control = list(indicator = 'inc',
                           control = list(adaptive = TRUE, adaptInterval = 100,
                                  scale = .5, log = FALSE, reflective = FALSE)))

mcmc <- buildMCMC(conf)
cmcmc <- compileNimble(mcmc, project = model)

cmcmc$run(5000)

smpDA <- as.matrix(cmcmc$mvSamples)

beta_enforce <- smpDA[ , 'beta_enforce']

par(mfrow = c(2,3))
ts.plot(smpDA[ , 'beta_benefit'])
ts.plot(smpDA[ , 'beta_age'])
ts.plot(smpDA[ , 'sigma_alpha'])
ts.plot(beta_enforce)
## beta_enforce when it is in the model
ts.plot(beta_enforce[beta_enforce != 0])  
hist(beta_enforce[beta_enforce != 0])

burnin <- 1000
smpDA <- smpDA[(burnin+1):nrow(smpDA), ]
## inclusion probability for enforcement
mean(smpDA[ , 'inc'])


