###############################################################
#   Simulation of the Wiener process on the interval [t0,T]
#                 Author : KACEF,M.A. PhD.
#                    ENSSEA, 2024
#                 
# X_t is a Wiener process on [t0,T] and starts from X0 if
#   1. X_t0=X0, i.e. the process starts from an initial value X0 at t0
#   2. The increments X_t(i)-X_t(i-1) are independent for i=1,...,n. 
#   3. For all t>s>=t0, we have X_t-X_s~N(C,sigma^2*(t-s)), where sigma>0 (volatility)
# Note: If sigma=1, X_t is called standard Brownian motion (SBM).
# 
################################################
#                  Arguments 
################################################
# X0 start value of the process, i.e. X(0)=X0
# The sigma	the volatility of this is a strictly positive numerical value; for example, 0.25 means a volatility of 25%.
# [t0,T] Time interval.
# N number of simulated paths
# n number of points used to construct a one path of X
# Z is a vector of n random variables distributed as N(0,1)
################################################

sbm=function(n,N,t0,T,X0,mu,sigma,Z){
  delta=(T-t0)/n
  X=matrix(nrow=N,ncol=n+1)
  X[,1]=X0
    for(i in 1:N){
    for(j in 2:(n+1)){
      X[i,j]=X[i,j-1]+sigma*sqrt(delta)*Z[i,j-1]
    }}  
  return(X)
}

#########################################
#           Plot paths function
########################################
plot_paths=function(X,t0,T)
{
  dt=(T-t0)/n
  t=seq(t0,T,by=dt) #length(t)=n+1
  plot(NA,type="n",xlim=c(t0,T),ylim=c(min(X),max(X)),ylab=expression(X[t]),xlab="Time",main=expression("Simulation of the sample paths"),cex.lab=1.2,cex.axis=1.2)
  for (i in 1:N)
  {
    lines(t,X[i,],col="#2297E6")  #length(X[i,])=n+1 for all i=1,...,N
  }
}

#########################################################
# Notes.
# -X_t(j)=X[,j] column vector of all N realisations
#  of the process at time t_j.
# -X_t(j)(i)=X[i,j] is the realisation of process X at time t_j for path i. 
# -The index j corresponds to the time position, i.e. j=t_{j-1}, for j=1,...,n+1,
#  So, for example, the instant at position j is equal to t0+(j-1)*delta 
#########################################################
# Simulation of the sample paths of SBM
#########################################################

n=10^3
N=10^2
t0=0
T=1
X0=0
sigma=1
Z=matrix(rnorm(N*n),nrow=N,ncol=n)
X=sbm(n,N,t0,T,X0,mu,sigma,Z)
plot_paths(X,t0,T)

# Mean function of X

dt=(T-t0)/n
t=seq(t0,T,by=dt)
lines(t,colMeans(X),col='red',lwd=2) # Empirical mean function 
abline(h=X0,lwd=2,lty=2)

