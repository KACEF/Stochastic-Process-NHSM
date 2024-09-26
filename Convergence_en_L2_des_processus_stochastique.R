###########################################
#  Converenge en moyenne quadratique 
#  Pont Brownien vs. Ornstein-Uhlenbeck 
###########################################

###########################################################
# Simulation du processus de Ornstein-Uhlenbeck (O-U) (1930)
#                Auteur. KACEF,MA, PhD.
#                 ENSSEA, mai 2023
##############################################################
#         dX(t)=theta(mu-X(t))dt+sigmadB_t
#           theta>0, mu in R, sigma>0
##############################################################


################################################
#          Arguments 
################################################
# X0 start value of the Arithmetic Brownian Motion, i.e. X(0)=X0
# # sigma	the volatility , this is a strictly positive numerical value;
# [t0,T] Time interval
# N number of simulated paths
# n number of points used to construct a one path of X
################################################
# Simulation of paths of U-O process in [t0,T]
#####################################################

uo=function(n,N,t0,T,X0,theta,mu,sigma,Z){
  delta=(T-t0)/n
  X=matrix(nrow=N,ncol=n+1)
  X[,1]=X0
  for(j in 1:N){
    for(i in 2:(n+1)){
      X[j,i]=X[j,i-1]+theta*(mu-X[j,i-1])*delta+sigma*sqrt(delta)*Z[j,i-1]
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
  plot(t,type="n",xlim=c(t0,T),ylim=c(min(X),max(X)),ylab=expression(X[t]),xlab="Time",main=expression("The sample paths of the U-O process "),cex.lab=1.2,cex.axis=1.2)
  for (j in 1:N)
  {
    lines(t,X[j,],col="#2297E6")
  }
}
########################################
# Example of paths simulation of U-O
########################################
n=10^3
N=50
t0=0
T=1
X0=0
mu=0 #moyenne à long terme 
sigma=1
theta<-2.75
Z=matrix(rnorm(N*n),nrow=N,ncol=n)
X=uo(n,N,t0,T,X0,theta,mu,sigma,Z)
par(mfrow=c(1,2))
plot_paths(X,t0,T)
abline(h=mu,lwd=2,col="green",lty=2)
###################################
# plot empirical mean of X
###################################
y<-vector()
for(i in 1:(n+1)){
  y[i]<-mean(X[,i])
}
t=seq(from=0,to=1,length.out=n+1)
lines(t,y,col="red",lwd=3,type='l',lty=2)
###########################################
#         Mean function of U-O process
#########################################
mean_UO<-function(t){X0*exp(-theta*t)+mu*(1-exp(-theta*t))}
# plot mean function of X
curve(mean_UO,from=0,to=1,n=10^3,col="black",add=TRUE,lwd=2)

#####################################################################
legend(x=0.6,y=3,legend=c("The mean function ","The empirical mean"), col =c("black","red"),lty=c(1,2),bty="n",lwd=c(3,3))
###################################################################

###################################
# plot empirical variance of X
###################################
y2<-vector()
for(i in 1:(n+1)){
  y2[i]<-var(X[,i])
}
t=seq(from=0,to=1,length.out=n+1)
lines(t,y2,col="yellow",lwd=3,type='l',lty=1)
###########################################
###########################################
#         Variance function of U-O process
#########################################
var_UO<-function(t){sigma^2/(2*theta)*(1-exp(-2*theta*t))}
# plot variance function of X
curve(var_UO,from=0,to=1,n=10^3,col="black",add=TRUE,lwd=2)

############################################
abline(h=sigma^2/theta*0.5,lwd=2,col="green",lty=2)
#####################################
#   Simulation d'un pont brownien 
#            Xt=Bt-t/T*BT
######################################
BB=function(n,N,t0,T,X0,Z){
  delta=(T-t0)/n
  X=matrix(nrow=N,ncol=n+1)
  X[,1]=X0
  for(j in 1:N){
    for(i in 2:(n+1)){
      X[j,i]=X[j,i-1]+delta+sqrt(delta)*Z[j,i-1]
    }}
  Y=matrix(nrow=N,ncol=n+1)
  Y[,1]=X0
  for(j in 1:N){
    for(i in 2:(n+1)){
      Y[j,i]=X[j,i]-i*delta/T*X[j,n+1]
    }}
return(Y)
}
#########################################
#           Plot paths function
########################################
plot_paths=function(X,t0,T)
{
  dt=(T-t0)/n
  t=seq(t0,T,by=dt)
  plot(t,type="n",xlim=c(t0,T),ylim=c(min(X),max(X)),ylab=expression(X[t]),xlab="Time",main=expression("The sample paths of Brownien Bridge process"),cex.lab=1.2,cex.axis=1.2)
  for (j in 1:N)
  {
    lines(t,X[j,],col="#2297E6")
  }
}
########################################
# Example of paths simulation of abm
########################################
n=10^3
N=50
t0=0
T=1
X0=0
Z=matrix(rnorm(N*n),nrow=N,ncol=n)
X=BB(n,N,t0,T,X0,Z)
plot_paths(X,t0,T)
################################################


###################################
# plot empirical mean of X
###################################
y<-vector()
for(i in 1:(n+1)){
  y[i]<-mean(X[,i])
}
t=seq(from=0,to=1,length.out=n+1)
lines(t,y,col="red",lwd=3,type='l',lty=2)
###########################################

###################################
# plot empirical variance of X
###################################
y2<-vector()
for(i in 1:(n+1)){
  y2[i]<-var(X[,i])
}
t=seq(from=0,to=1,length.out=n+1)
lines(t,y2,col="yellow",lwd=3,type='l',lty=1)
###########################################
abline(h=0,lwd=2,col="black",lty=1)

###########################################
#         Variance function of BB process
#########################################
var_BB<-function(t){t-t^2/T}
# plot variance function of X
curve(var_BB,from=0,to=1,n=10^3,col="black",add=TRUE,lwd=2)

