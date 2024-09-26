#################################
#         Random Walk
#    Author. KACEF,MA, PhD
#           ENSSEA
#################################
#
#    S(t)=sum Xi, S0=0
#     X=+1 with p
#     X=-1 with 1-p  
#     p in ]0,1[
#################################

X_i<-function(p){
  u<-runif(1)
  if (u<p) return(+1)
  else return(-1)}

# Simulation of one path of S

  S<-function(p,t0,T,S0){
  n<-T-t0
  X<-vector()
  X<-replicate(n,X_i(p))
  S<-vector()
  S[1]<-S0
  S<-c(S0,cumsum(X))
  return(S)
}

  #####################################
  #             INPUTS
  #####################################
  # S0 start value of the AR(1), i.e. S(0)=S0.
  # p is the probability of taking a step forward, i.e. X=+1
  # [t0,T] Time interval.
  # N number of simulated paths
  ###################################### 

N<-50
p<-0.53
S0<-0
t0<-0
T<-10^2

##### Simulation ####################

DATA<-replicate(N,S(p,t0,T,S0))

# plot sample paths of S in [t0,T]

t<-seq(from=t0, to=T, by=1)
plot(t,type="n",ylab=expression(S[t]),xlab="t",ylim=c(min(DATA),max(DATA)),main="Sample paths for a random walk",xlim=c(t0,T),cex.lab=1.2,cex.axis=1.2)
for( i in 1:N)
{lines(t,DATA[,i],col="#2297E6",lwd=2)}


######### plot mean of S ##############

y<-vector()
for(i in 1:length(t)){
  y[i]<-mean(DATA[i,])
}

lines(t,y,col="red",lwd=3,type='l')
legend(x=0,y=max(DATA),legend =expression("Mean of process "* S[t]), col ="red",lty=1,bty="n",lwd=3)

##########################################
#      True mean of process S
##########################################

mean_rw<-function(t) {(2*p-1)*t}
y2<-vector()
y2<-sapply(t,mean_rw)
lines(t,y2,col="black",lwd=3,type='l')
legend(x=0,y=max(DATA)-3,legend =expression(E(S[t])), col ="black",lty=1,bty="n",lwd=3)
##########################################