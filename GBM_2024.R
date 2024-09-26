################################################
#   Simulation of Geometric Brownian motion paths
#     SDE: dSt=mu*S_t*dt+sigma*S_t*dBt,
# where, S0>0, mu in R and sigma>0 (volatility)
# The unique solution of this SDE is obtained with Ito's formula as
# S_t=S0*exp((mu-0.5*sigma^2)*t+sigma*B_t)
# # Remark. 
# If X is an arithmetic Brownian motion with drift mu-0.5*sigma^2 and volatility sigma, then 
# S_t=S0*exp(X_t)
# ##
# Author : KACEF,M.A.
# ENSSEA, 2024
################################################
#                  Arguments 
################################################
# S0 start value of the process, i.e. S(0)=S0>0
# The sigma	the volatility of this is a strictly positive numerical value; for example, 0.25 means a volatility of 25%.
# [t0,T] Time interval.
# N number of simulated paths
# n number of points used to construct a one path of S
# Z is a vector of n random variables distributed as N(0,1)
################################################

gbm=function(n,N,t0,T,S0,mu,sigma,Z){
  delta=(T-t0)/n
  X=matrix(nrow=N,ncol=n+1)
  X[,1]=S0
  for(i in 1:N){
    for(j in 2:(n+1)){
      X[i,j]=X[i,j-1]+mu*X[i,j-1]*delta+sigma*X[i,j-1]*sqrt(delta)*Z[i,j-1]
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
# Example of paths simulation of GBM
########################################
########################################

n=10^3
N=50
t0=0
T=1
dt=(T-t0)/n
S0=1
mu=1
sigma=0.25
Z=matrix(rnorm(N*n),nrow=N,ncol=n)
X=gbm(n,N,t0,T,S0,mu,sigma,Z)
plot_paths(X,t0,T)

#########################################
# Mean function of X
# #######################################
# Empirical mean function 
t=seq(t0,T,by=dt)
lines(t,colMeans(X),col='red',lwd=2) 
# Theorical mean function
m=function(t){S0*exp(mu*t)}
curve(m,from=t0,to=T,lwd=2,n=10^3,add=TRUE,lty=2)



