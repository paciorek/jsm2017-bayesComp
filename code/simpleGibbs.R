### Data on body temperature for 12 people
y = c(98.0, 97.5, 98.6, 98.8, 98.0, 98.5, 98.6, 99.4, 98.4, 98.7, 98.6, 97.6)
n = length(y)

### Prior parameters
mu    = 98
tau2  = 100
aa    = 2
cc    = 1

### Parameter initialization
theta  = theta0  = 98
sigma2 = sigma20 = 1

### Vectors to store
B = 11000
theta.out = rep(0, B)
sigma2.out = rep(0, B)

### Gibbs sampler
for(b in 1:B){
  theta  = rnorm(1, (sum(y)/sigma2 + mu/tau2)/(n/sigma2 + 1/tau2), sqrt(1/(n/sigma2 + 1/tau2)))  #Remember R parameterizes the univariate normal in terms of the standard deviation
  sigma2 = 1/rgamma(1, aa+n/2, rate=cc+sum((y-theta)^2)/2)  #To avoid mistakes, I always tell R to use the rate rather than the scale ...

  ### Decoupling the current value from the parameters from where it is stored is very helpful for debugging, I strongly suggest you always do this!
  theta.out[b]  = theta
  sigma2.out[b] = sigma2
  if(b/50==floor(b/50)){
    print(b)
  }
}

### Trace plots
quartz()
plot(theta.out, type="l")
quartz()
plot(sigma2.out, type="l")

### Posterior distributions
quartz()
plot(theta.out, sigma2.out)
quartz()
hist(theta.out)
abline(v=theta0, col="red")
quartz()
hist(sigma2.out)
abline(v=sigma20, col="red")

theta.out.fin  = theta.out[-seq(1,1000)]
sigma2.out.fin = sigma2.out[-seq(1,1000)]
save(y, theta.out.fin, sigma2.out.fin, file="simpleGibbsoutput.Rdata")
