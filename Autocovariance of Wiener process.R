
##########################
#  Autocovariance of Wiener process with
#  Author. KACEF,MA. PhD.
#  ENSSEA, 2023
###########################

################################################
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
  plot(t,type="n",xlim=c(t0,T),ylim=c(min(X),max(X)),ylab=expression(X[t]),xlab="Time",main=expression("Paths of a Wiener process with sigma=0.5"),cex.lab=1.2,cex.axis=1.2)
  for (j in 1:N)
  {
    lines(t,X[j,],col="#2297E6")
  }
  
}
########################################
# Example of paths simulation of abm
########################################
n=10^2
N=10^3
t0=0
T=1
X0=1
mu=0
sigma=0.5
Z=matrix(rnorm(N*n),nrow=N,ncol=n)
X=abm(n,N,t0,T,X0,mu,sigma,Z)
par(mfrow=c(1,2))
plot_paths(X,t0,T)

##########################################
# 3 rd Qu.
y<-vector()
for(j in 1:(n+1)) {y[j]<-mean(X[,j])}
dt=(T-t0)/n
t=seq(t0,T,by=dt)
lines(t,y,lwd=3,col="red")
###########################################
legend(x=0,y=2.9,legend =expression("Mean of process "* X[t]), col ="red",lty=1,bty="n",lwd=3)
##########################################
# True mean of process S

C_x<-function(t) {sigma^2*t}
y2<-vector()
y2<-sapply(t,C_x)
plot(t,y2,col="black",lwd=3,type='l',xlim=c(t0,T),ylim=c(0,sigma^2),ylab="The autocovariance function",xlab="Time",main=expression(""),cex.lab=1.2,cex.axis=1.2)
legend(x=0,y=,0.25,legend =expression(Cov(X[t],X[t+s])), col ="black",lty=1,bty="n",lwd=3)
legend(x=0,y=,0.235,legend =expression("Approx."), col ="springgreen3",lty=1,bty="n",lwd=3)


######################

covat<-vector()
for(i in 1:n){
  covat[i]<-cov(X[,i],X[,i+1])
}
y3<-vector()
y3<-covat
t3<-t[-(n+1)]
lines(t3,y3,col="springgreen3",lwd=3,type='l')
