RW_sampler_nonzero_indicator <- nimbleFunction(
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
        regular_RW_sampler <- sampler_RW(model, mvSaved, target = target, control = control$control)
        indicatorNode <- control$indicator
    },
    run = function() {
        if(model[[indicatorNode]] == 1) regular_RW_sampler$run()
    },
    methods = list(
        reset = function() {regular_RW_sampler$reset()}
    ))

RJindicatorSampler <- nimbleFunction(
    contains = sampler_BASE,
    setup = function( model, mvSaved, target, control ) {
        ## target should be the name of the indicator node, 'z2' above
        ## control should have an element called coef for the name of the corresponding coefficient ('beta2' above.  This could potentially be determined from the model structure, but for now I won't try that.)
        coefNode <- control$coef
        scale <- control$scale
        calcNodes <- model$getDependencies(c(coefNode, target))
        ## coefNode not in reduced model so its prior not calculated
        calcNodesReduced <- model$getDependencies(target)
    },
    run = function( ) {
        currentIndicator <- model[[target]]
        if(currentIndicator == 1) {
            ## propose removing it
            currentLogProb <- model$getLogProb(calcNodes)
            currentCoef <- model[[coefNode]]
            ## reverse jumping density (density for proposed auxiliary variable (i.e., current coefficient)
            logProbReverseProposal <- dnorm(currentCoef, 0, sd = scale, log = TRUE)
            model[[target]] <<- 0
            model[[coefNode]] <<- 0
            model$calculate(calcNodes)
            ## avoid including prior for coef not in model
            log_accept_prob <- model$getLogProb(calcNodesReduced) - currentLogProb + logProbReverseProposal
        } else {
            ## propose adding it
            currentLogProb <- model$getLogProb(calcNodesReduced)
            proposalCoef <- rnorm(1, 0, sd = scale)
            model[[target]] <<- 1
            model[[coefNode]] <<- proposalCoef
            ## jumping density (density for current auxiliary variable (i.e., proposed coefficient)
            logProbForwardProposal <- dnorm(proposalCoef, 0, sd = scale, log = TRUE)
            proposalLogProb <- model$calculate(calcNodes)
            log_accept_prob <- proposalLogProb - currentLogProb - logProbForwardProposal
        }
        accept <- decide(log_accept_prob)
        if(accept) {
            copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
        } else {
            copy(from = mvSaved, to = model, row = 1, nodes = calcNodes, logProb = TRUE)
        }
    },
    methods = list(reset = function() {
    })
)

