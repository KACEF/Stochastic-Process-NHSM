##################################################
#    The Ornstein-Uhlenbeck process 
#             KACEF, MA. PhD.
#                MAY 24'
##################################################

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
  plot(NA,type="n",xlim=c(t0,T),ylim=c(min(X),max(X)),ylab=expression(X[t]),xlab="Time",main="",cex.lab=1.2,cex.axis=1.2)
  title(main = list("Simulation of the sample paths of O-U process", cex =0.9,
                    col = "black", font = 1))
  for (i in 1:N)
  {
    lines(t,X[i,],col="#2297E6")  #length(X[i,])=n+1 for all i=1,...,N
  }
}

#########################################################
#           Simulation of Ornstein-Uhlenbeck 
#          process (the mean-reverting process)
#########################################################

n=10^3
N=500
t0=0
T=1
X0=50
mu=5        # The long-term mean
sigma=8     # Volatility of the process.
theta<-2    # The speed of mean-reverting.
alpha=theta*mu
beta=-theta
gamma=0
Z=matrix(rnorm(N*n),nrow=N,ncol=n)
X=CKLS(n,N,t0,T,X0,alpha,beta,gamma,sigma,Z)
par(mfrow=c(1,2))
plot_paths(X,t0,T)
abline(h=mu,lwd=2,lty=1,col='green')

# Analysis of mean 

dt=(T-t0)/n
t=seq(t0,T,by=dt)
lines(t,colMeans(X),col='red',lwd=3,lty=3) # Empirical mean function 
m=function(t){(X0-mu)*exp(-theta*t)+mu} # Theorical mean function
curve(m,from=t0,to=T,lwd=2,n=10^3,add=TRUE)
legend(x=0.4,y=50,legend=c("Empirical mean","Long-term mean","Theorical mean function"), col =c("red","green",'black'),cex=0.5,lty=c(3,1,1),bty="n",lwd=c(3,2,2))

# Analysis of variance

plot(t,colMeans(X^2)-colMeans(X)^2,col='red',lwd=2,xlab="Time",type='l',ylab=expression(Var*(X[t])),main="") # Empirical variance
title(main = list("Variance of O-U process", cex =0.9,
                  col = "black", font = 1))
var=function(t){0.5*sigma^2/theta*(1-exp(-2*theta*t))} # Theorical variance
curve(var,from=t0,to=T,lwd=2,n=10^3,add=TRUE)
legend(x=0.4,y=10,cex=0.5,legend=c("Empirical variance","Long-term variance","Theorical variance"), col =c("red","green","black"),lty=c(1,1,1),bty="n",lwd=c(2,2,2))
abline(h=sigma^2*0.5/theta,lwd=2,lty=1,col='green')

