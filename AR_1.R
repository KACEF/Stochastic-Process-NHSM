################################
#  Simulation of AR(1)
#  Author. KACEF, MA.PhD
#  ENSSEA, 2023
################################


  ar_1=function(n,N,t0,T,X0,mu,sigma){
    delta=1
    X=matrix(nrow=N,ncol=n+1)
    X[,1]=X0
    for(j in 1:N){
    for(i in 2:(n+1)){
      Z<-vector()
      Z<-rnorm(n=i-1,mean=mu,sd=sigma)
      X[j,i]=X0+sum(Z)
      }}  
    return(X)
  }
  
  #########################################
  #           Plot paths function
  ########################################
  plot_paths=function(X,t0,T)
  {
    dt=1
    t=seq(t0,T,by=dt)
    plot(t,type="n",xlim=c(t0,T),ylim=c(min(X),max(X)),ylab=expression(X[t]),xlab="Time",main=expression("Path of an AR(1) with mu=1 and sigma=0.48"),cex.lab=1.2,cex.axis=1.2)
    for (j in 1:N)
    {
      lines(t,X[j,],col="#2297E6",type='l',lwd=2)
    }
    
  }
  
  #####################################
  #             INPUTS
  #####################################
  # X0 start value of the AR(1), i.e. X(0)=X0.
  # mu the drift parameter. 
  # sigma	the volatility , this is a strictly positive numerical value; for example, 0.25 means a volatility of 25%.
  # [t0,T] Time interval.
  # N number of simulated paths
  # n number of points used to construct a one path of AR(1)
  ################################################ 
  
  n<-100
  N<-2
  mu<-2
  sigma<-025
  t0<-0
  T<-n
  X0<-1
  
  ########## Simulation ####################  
  
  X<-ar_1(n,N,t0,T,X0,mu,sigma)
  plot_paths(X,t0,T)
  
  ######################################### 