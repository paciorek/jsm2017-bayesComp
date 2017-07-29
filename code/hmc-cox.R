library(rstan)

## read in data and initial values
data = new.env()
source('leuk.data.R', local = data)
data <- as.list(data)
init = new.env()
source('leuk.init.R', local = init)
init <- as.list(init)

## fit model -- initial delay is generation and compilation of C++
fit <- stan('leuk.stan', chains = 1, data = data, init = list(init))

## subsetting is because we only run one chain
smp <- extract(fit, permuted = FALSE, inc_warmup = TRUE)[ , 1, ]

ts.plot(smp[ , 'beta'])
coda::effectiveSize(smp[100:2000, 'beta'])
          
hazCols <- paste0('dL0[', 1:data$NT, ']')
SplaceboCols <- paste0('S_placebo[', 1:data$NT, ']')
StreatCols <- paste0('S_treat[', 1:data$NT, ']')

## plotting of estimated survival functions
mn <- colMeans(smp[100:2000,])
tmp <- c(0, data$t)
len <- length(tmp)
midPoints <- (tmp[1:(len-1)] + tmp[2:len]) / 2
plot(midPoints, c(1, mn[SplaceboCols]), ylim = c(0, 1), col = 'darkgrey',
     xlab = 'time (days)', ylab = 'survival probability')
points(midPoints, c(1, mn[StreatCols]))
legend('topright', col = c('darkgrey', 'black'),
       legend = c('placebo', 'treatment'))
