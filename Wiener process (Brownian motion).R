###########################################
#   Simulation of arithmetic Brownian motion paths
#        SDE:  dXt=mu*dt+sigma*dWt
#   mu in R (drift) and  Sigma>0 (volatility)
#            Author: KACEF,M.A.
#              ENSSEA, 2022
################################################

################################################
#          Arguments 
################################################
# X0 start value of the Arithmetic Brownian Motion, i.e. X(0)=X0
# mu the drift parameter 
# sigma	the volatility of the underlying asset, this is a strictly positive numerical value; for example, 0.25 means a volatility of 25%.
# [t0,T] Time interval, in option pricing we consider that t0 is the purchase date of the contract and T is its maturity.
# N number of simulated paths
# n number of points used to construct a one path of X
################################################

abm=function(n,N,t0,T,X0,mu,sigma,Z){
  delta=(T-t0)/n
  X=matrix(nrow=N,ncol=n+1)
  X[,1]=X0
  for(j in 1:N){
    for(i in 2:(n+1)){
      X[j,i]=X[j,i-1]+mu*delta+sigma*sqrt(delta)*Z[j,i-1]
    }}  
  return(X)
}

#########################################
#           Plot paths function
########################################
plot_paths=function(X,t0,T)
{
  dt=(T-t0)/n
  t=seq(t0,T,by=dt)
  plot(t,type="n",xlim=c(t0,T),ylim=c(min(X),max(X)),ylab=expression(X[t]),xlab="Time",main=expression("The sample paths of the Wiener process with sigma=0.25 and X0=1"),cex.lab=1.2,cex.axis=1.2)
  for (j in 1:N)
  {
    lines(t,X[j,],col="#2297E6")
  }
}
########################################
# Example of paths simulation of abm
########################################
n=10^3
N=100
t0=0
T=1
X0=0
mu=0
sigma=1
Z=matrix(rnorm(N*n),nrow=N,ncol=n)
X=abm(n,N,t0,T,X0,mu,sigma,Z)
plot_paths(X,t0,T)

# plot mean of S
y<-vector()
for(i in 1:(n+1)){
  y[i]<-mean(X[,i])
}
t=seq(from=0,to=1,length.out=n+1)
lines(t,y,col="red",lwd=3,type='l')
legend(x=0,y=1.6,legend =expression("Mean of process "* X[t]), col ="red",lty=1,bty="n",lwd=3)
##########################################




