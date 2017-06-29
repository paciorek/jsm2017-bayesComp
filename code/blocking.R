library(mvtnorm)

rrr  = 11000
burn = 1000
sampler1 = function(rrr, burn){
  theta.out = array(NA, dim=c(rrr-burn, 3))
  theta = c(2, 1, -1)
  for(s in 1:rrr){
    theta[1] = rgamma(1, 3, rate=0.5*(1 + theta[2]^2 + theta[3]^2 - 2*theta[2]*theta[3]/sqrt(2)))
    theta[2] = rnorm(1, theta[3]/sqrt(2), 1/sqrt(theta[1]))
    theta[3] = rnorm(1, theta[2]/sqrt(2), 1/sqrt(theta[1]))

    if(s > burn){
      theta.out[s-burn,] = theta
    }
  }
  return(theta.out)
}

sampler2 = function(rrr, burn){
  theta.out = array(NA, dim=c(rrr-burn, 3))
  theta = c(2, 1, -1)
  AA = matrix(c(2, sqrt(2), sqrt(2), 2), 2, 2)
  for(s in 1:rrr){
    theta[1] = rgamma(1, 3, rate=0.5*(1 + theta[2]^2 + theta[3]^2 - 2*theta[2]*theta[3]/sqrt(2)))
    theta[2:3] = rmvnorm(1, rep(0,2), AA/theta[1])
    if(s > burn){
      theta.out[s-burn,] = theta
    }
  }
  return(theta.out)
}


#################################################################
res1 = sampler1(rrr, burn)
quartz()
par(mfrow=c(2,2))
par(mar=c(4,4,1,1)+0.2)
plot(res1[,1], type="l", ylab=expression(theta[1]))
abline(h = 4, col="grey")
par(mar=c(4,4,1,1)+0.2)
acf(res1[,1])
#
par(mar=c(4,4,1,1)+0.2)
plot(res1[,2], type="l", ylab=expression(theta[2]))
abline(h = 0, col="grey")
par(mar=c(4,4,1,1)+0.2)
acf(res1[,2])
dev.print(dev=pdf, file="blocking1.pdf")



res2 = sampler2(rrr, burn)
quartz()
par(mfrow=c(2,2))
par(mar=c(4,4,1,1)+0.2)
plot(res2[,1], type="l", ylab=expression(theta[1]))
abline(h = 4, col="grey")
par(mar=c(4,4,1,1)+0.2)
acf(res2[,1])
#
par(mar=c(4,4,1,1)+0.2)
plot(res2[,2], type="l", ylab=expression(theta[2]))
abline(h = 0, col="grey")
par(mar=c(4,4,1,1)+0.2)
acf(res2[,2])
dev.print(dev=pdf, file="blocking2.pdf")







