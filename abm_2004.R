################################################
#   Simulation of arithmetic Brownian motion paths
#        SDE:  dXt=mu*dt+sigma*dBt
#   mu in R (drift) and  Sigma>0 (volatility)
#            Author: KACEF,M.A.
#              ENSSEA, 2024
################################################

################################################
#                  Arguments 
################################################
# X0 start value of the process, i.e. X(0)=X0
# The sigma	the volatility of this is a strictly positive numerical value; for example, 0.25 means a volatility of 25%.
# [t0,T] Time interval.
#  N number of simulated paths
# n number of points used to construct a one path of X
# Z is a vector of n random variables distributed as N(0,1)
################################################

abm=function(n,N,t0,T,X0,mu,sigma,Z){
  delta=(T-t0)/n
  X=matrix(nrow=N,ncol=n+1)
  X[,1]=X0
  for(i in 1:N){
  for(j in 2:(n+1)){
    X[i,j]=X[i,j-1]+mu*delta+sigma*sqrt(delta)*Z[i,j-1]
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
  plot(t,type="n",xlim=c(t0,T),ylim=c(min(X),max(X)),ylab=expression(X[t]),xlab="Time",main=expression("Simulation of the sample paths"),cex.lab=1.2,cex.axis=1.2)
  for (j in 1:N)
  {
    lines(t,X[j,],col="#2297E6")
  }
}
########################################
# Example of paths simulation of abm
########################################
n=10^3
N=10^2
t0=0
T=1
dt=(T-t0)/n
X0=2
mu=3.25
sigma=1
Z=matrix(rnorm(N*n),nrow=N,ncol=n)
X=abm(n,N,t0,T,X0,mu,sigma,Z)
plot_paths(X,t0,T)


#########################################
# Mean function of X
# #######################################
# Empirical mean function 
t=seq(t0,T,by=dt)
lines(t,colMeans(X),col='red',lwd=2) 
# Theorical mean function
m=function(t){mu*t+X0}
curve(m,from=t0,to=T,lwd=2,n=10^3,add=TRUE,lty=2)


