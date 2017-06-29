lambda.true = 2
n.s = c(10, 50)
for(i in 1:length(n.s)){
  n = n.s[i]
  y = rep(lambda.true, n)
  xx = seq(0,6,length=100)
  quartz()
  plot(xx, dgamma(xx, a + sum(y), rate=b + n),type="l", xlab=expression(lambda), ylab="Density", lwd=2, main=paste("n =",n))
  lines(xx, dnorm(xx, (a - 1 + sum(y))/(b + n), sqrt((a - 1 + sum(y))/((b + n)^2))), col="red", lwd=2, lty=2)
  lines(xx, dnorm(xx, mean(y), sqrt(mean(y)/n)), col="blue", lwd=2, lty=3)
  legend(4, 0.95*max(dgamma(xx, a + sum(y), rate=b + n)), c("Exact", "Like + Prior", "Like only"), lty=c(1,2,3), lwd=2, col=c("black", "red", "blue"), bty="n")
  dev.print(dev=pdf, file=paste("laplaceapprox_poissongamma_n",n,".pdf", sep=""))
}


## Posterior means
lambda.true = 2
n.s = c(10, 50, 500)
a.s = c(1, 2)
b.s = c(1, 1/5)

for(i in 1:length(n.s)){
  n = n.s[i]
  y = rep(lambda.true, n)
  for(j in 1:length(a.s)){
    a = a.s[j]
    b = b.s[j]
    
    print(paste("Posterior means", "n =", n,"    a =", a,"    b =", b))
    print(paste("Exact ",(a + sum(y))/(b + n)))
    print(paste("Like + prior ",(a - 1 + sum(y))/(b + n)))
    print(paste("Like ",sum(y)/n))

  }
}



## Tail areas
n.s  = c(10, 50)
ll.s = c(1.5, 1.8, 2.5)
ul.s = c(2.8, 2.3, Inf)
a    = 2
b    = 1/5


for(i in 1:length(n.s)){   #
  n = n.s[i]
  y = rep(lambda.true, n)
  for(j in 1:length(ll.s)){   #
    ll = ll.s[j]
    ul = ul.s[j]
    
    print(paste("Posterior areas", "n =", n,"    ll =", ll,"    ul =", ul))
    print(paste("Exact ", pgamma(ul, a + sum(y), rate=b + n) - pgamma(ll, a + sum(y), rate=b + n) ))
    
    
    
    print(paste("Like + prior ", pnorm(ul, (a - 1 + sum(y))/(b + n), sqrt((a - 1 + sum(y))/((b + n)^2))) - pnorm(ll, (a - 1 + sum(y))/(b + n), sqrt((a - 1 + sum(y))/((b + n)^2))) ))
    print(paste("Like ",  pnorm(ul, mean(y), sqrt(mean(y)/n)) - pnorm(ll, mean(y), sqrt(mean(y)/n)) ))
  }
}



