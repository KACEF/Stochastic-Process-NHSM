#######################################
#     Simulation of process 
#    Xt=epsilon1+epsilon2*t^2
######################################
n<-10^2
N=10^2
lambda=0.75
t<-seq(from=0,to=1,length.out=n)
x_t<-function(t) {1+rexp(1,rate=lambda)*t^2}
DATA<-replicate(N,x_t(t))
plot(t,type="n",ylab=expression(X[t]),xlab="t",ylim=c(min(DATA),max(DATA)),main="Sample paths for process X (Lambda=0.75, N=100)",xlim=c(t0,T),cex.lab=1.2,cex.axis=1.2)
for( i in 1:N)
{lines(t,DATA[,i],col="#2297E6",lwd=2)}


# True mean of process X

mean_rw<-function(t) {1+1/lambda*t^2}
y2<-vector()
y2<-sapply(t,mean_rw)
lines(t,y2,col="red",lwd=3,type='l')
legend(x=0,y=6,legend =expression(E(X[t])), col ="red",lty=1,bty="n",lwd=3)
##########################################