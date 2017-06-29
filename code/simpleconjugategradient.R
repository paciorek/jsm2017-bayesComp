### Data on body temperature for 12 people
y = c(98.0, 97.5, 98.6, 98.8, 98.0, 98.5, 98.6, 99.4, 98.4, 98.7, 98.6, 97.6)
n = length(y)

### Prior parameters
mu    = 95
tau2  = 1
aa    = 2
cc    = 1

### Parameter initialization
theta  = theta.o  = 98
sigma2 = sigma2.o = 1

### Parameters controlling algorithm stopping.
maxiter = 200
epsilon = 10e-6
s       = 1

### Conjugate gradient algorithm
while(s<=maxiter){
  theta  = (sum(y)/sigma2 + mu/tau2)/(n/sigma2 + 1/tau2)
  sigma2 = (cc+sum((y-theta)^2)/2)/(aa+1+n/2)
  if(max(abs(theta-theta.o)/theta.o, abs(sigma2.o-sigma2)/sigma2.o) < epsilon){
    break
  }
  s = s+1
  theta.o  = theta
  sigma2.o = sigma2
  print(paste(s, theta, sigma2))
}


print(paste("The MLE is: ", mean(y), (n-1)*var(y)/n))

