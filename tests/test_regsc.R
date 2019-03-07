library(regsc)

n=120
p=2
q=3
sigma=0.5

weight=rep(1,n-1)
XTol=1e-6
maxIter=1000

x=rnorm(n*p)
dim(x)<-c(n,p)
z=rnorm(n*q)
dim(z)<-c(n,q)

beta0=c(rep(1,n/2),rep(0.5,n/2))
beta0=rep(beta0,p)
dim(beta0)<-c(n,p)
gamma0=rep(1,q)
dim(gamma0)<-c(q,1)

y = rowSums(x*beta0) + z %*% gamma0 + sigma*rnorm(n)

data = as.data.frame(cbind(y,x,z))
colnames(data) <- c("y","x1","x2","z1","z2","z3")
  
#res=regsc(y~0+x1+x2|z1+z2+z3,data,method="rot")
res=regsc(y~0+x1+x2|0+z1+z2+z3,data)
summary(res)
#regime=post$regime
#print(regime)
#m=length(regime)-2
#alpha=post$alpha[1:(p*(m+1))]
#dim(alpha)<-c(p,m+1)
#alpha<-t(alpha)
#print(alpha)


