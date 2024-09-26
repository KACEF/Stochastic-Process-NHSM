##################################################################################
#            Simulation of the CKLS process (see [1])
#         SDE: dXt=(alpha+beta*Xt)*dt+sigma*Xt^gamma*dBt,
#     where, sigma>0,gamma>=0,alpha,beta in R and Bt is an SBM.
# Note.
# The Xt process is non-negative when alpha>0, beta<=0, gamma>0.5 and X0>0.
# References. 
# [1] CHAN, K. C., KAROLYI, G. A., LONGSTAFF, F. A., & SANDERS, A. B. (1992).
#  An Empirical Comparison of Alternative Models of the Short-Term Interest Rate.
#   The Journal of Finance, 47(3), 1209–1227.
#################################################################################
################################################
#                  Arguments 
################################################
# X0 start value of the process, i.e. X(0)=X0
# The sigma	the volatility of this is a strictly positive numerical value; for example, 0.25 means a volatility of 25%.
# [t0,T] Time interval, in option pricing we consider that t0 is the purchase date of the contract and T is its maturity.
# N number of simulated paths
# n number of points used to construct a one path of X
# Z is a vector of n random variables distributed as N(0,1)
################################################

CKLS=function(n,N,t0,T,X0,alpha,beta,gamma,sigma,Z){
  delta=(T-t0)/n
  X=matrix(nrow=N,ncol=n+1)
  X[,1]=X0
  for(i in 1:N){
    for(j in 2:(n+1)){
      X[i,j]=X[i,j-1]+(alpha+beta*X[i,j-1])*delta+sigma*(X[i,j-1])^gamma*sqrt(delta)*Z[i,j-1]
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
  plot(NA,type="n",xlim=c(t0,T),ylim=c(min(X),max(X)),ylab=expression(X[t]),xlab="Time",main=expression("Simulation of the sample paths of CKLS process"),cex.lab=1.2,cex.axis=1.2)
  for (i in 1:N)
  {
    lines(t,X[i,],col="#2297E6")  #length(X[i,])=n+1 for all i=1,...,N
  }
}

#########################################################
# Simulation of the sample paths of CKLS process
#########################################################

n=10^3
N=10^2
t0=0
T=1
X0=1.5
alpha=0.5
beta=-1.25
sigma=1
gamma=0.7
Z=matrix(rnorm(N*n),nrow=N,ncol=n)
X=CKLS(n,N,t0,T,X0,alpha,beta,gamma,sigma,Z)
plot_paths(X,t0,T)

# Mean function of X
dt=(T-t0)/n
t=seq(t0,T,by=dt)
lines(t,colMeans(X),col='red',lwd=2) # Empirical mean function 


#########################################################
# Example 1. Simulation of SBM
#########################################################

n=10^3
N=10^2
t0=0
T=1
X0=1.5
alpha=0
beta=0
sigma=1
gamma=0
Z=matrix(rnorm(N*n),nrow=N,ncol=n)
X=CKLS(n,N,t0,T,X0,alpha,beta,gamma,sigma,Z)
plot_paths(X,t0,T)

# Mean function of X
dt=(T-t0)/n
t=seq(t0,T,by=dt)
lines(t,colMeans(X),col='red',lwd=2) # Empirical mean function 

#########################################################
# Example 2. Simulation of ABM
#########################################################

n=10^3
N=10^2
t0=0
T=1
X0=1.5
alpha=2
beta=0
sigma=1.25
gamma=0
Z=matrix(rnorm(N*n),nrow=N,ncol=n)
X=CKLS(n,N,t0,T,X0,alpha,beta,gamma,sigma,Z)
plot_paths(X,t0,T)

# Mean function of X
dt=(T-t0)/n
t=seq(t0,T,by=dt)
lines(t,colMeans(X),col='red',lwd=2) # Empirical mean function 
# Theorical mean function
m=function(t){alpha*t+X0}
curve(m,from=t0,to=T,lwd=2,n=10^3,add=TRUE,lty=2)

#########################################################
# Example 3. Simulation of GBM
#########################################################

n=10^3
N=10^2
t0=0
T=1
X0=1.5
alpha=0
beta=1
sigma=0.25
gamma=1
Z=matrix(rnorm(N*n),nrow=N,ncol=n)
X=CKLS(n,N,t0,T,X0,alpha,beta,gamma,sigma,Z)
plot_paths(X,t0,T)

# Mean function of X
dt=(T-t0)/n
t=seq(t0,T,by=dt)
lines(t,colMeans(X),col='red',lwd=2) # Empirical mean function 
# Theorical mean function
m=function(t){X0*exp(beta*t)}
curve(m,from=t0,to=T,lwd=2,n=10^3,add=TRUE,lty=2)
legend(x=0,y=8,legend =expression("Mean of process "* X[t]), col ="red",lty=1,bty="n",lwd=2)

#########################################################
#           Simulation of Ornstein-Uhlenbeck 
#          process (the mean-reverting process)
#########################################################

n=10^3
N=10^2
t0=0
T=1
X0=1.5
mu=3         # The long-term average
sigma=0.5325 # Volatility of the process. Strictly positive. 
theta<-2.75  # The speed of mean-reverting. Strictly positive
alpha=theta*mu
beta=-theta
gamma=0
Z=matrix(rnorm(N*n),nrow=N,ncol=n)
X=CKLS(n,N,t0,T,X0,alpha,beta,gamma,sigma,Z)
plot_paths(X,t0,T)
abline(h=mu,lwd=2,lty=2)
# Mean function of X
dt=(T-t0)/n
t=seq(t0,T,by=dt)
lines(t,colMeans(X),col='red',lwd=2) # Empirical mean function 
legend(x=0.4,y=2,legend=c("Mean of process","Long-term mean"), col =c("red","black"),lty=c(1,2),bty="n",lwd=c(3,3))
###################################################################

