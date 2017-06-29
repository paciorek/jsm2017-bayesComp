###Importance sampling for the posterior mean of theta (prior is used as the importance distribution)
#  x_i | \theta ~ N(\theta, 1)
#  \theta ~ Cauchy (0, 1)
y = scan(file="impsampex_data.txt")

logposterior = function (y, theta){
  z = sum(dnorm(y, theta, 1, log=T)) + dcauchy(theta, log=T)
  return(z)
}

B        = 10000
theta.s1 = rcauchy(B, 0, 1)
lweights = rep(0, B)

for(b in 1:B){
  lweights[b] = logposterior(y, theta.s1[b]) - dcauchy(theta.s1[b], 0, 1, log=T)  #Note that the weights are computed on the log scale!!!!
}

weights = exp(lweights - max(lweights))
weights = weights/sum(weights)
quartz()
hist(weights)
quartz()
plot(theta.s1,weights, xlim=c(-50,50))
print(sum(weights*theta.s1))
print(sum(weights)^2/sum(weights^2))


