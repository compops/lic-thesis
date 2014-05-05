########################################################################
#
# Code to reproduce examples from the Licentiate's thesis:
# Sequential Monte Carlo for inference in nonlinear state space models
#
# written by
# Johan Dahlin ( johan.dahlin (at) liu.se
#
# Available from:
# http://users.isy.liu.se/en/rt/johda87/
#
# Copyright (c) 2014 Johan Dahlin
# 
# Run this code together with the corresponding R-file to reproduce:
#
# The IS algorithm for parameter inference
# Example 3.1 in Section 3.2 
#
########################################################################

th = matrix(0, 2000, 3)
ll = matrix(0, 2000, 1)

for ( ii in 0:7 ) {
  d <- read.table(paste(paste("ch3-example-importancesampling-hwsv",ii,sep=""),".csv",sep=""),header=T,sep=",",stringsAsFactors=F)[,-1]  
  th[ (ii*250+1):((ii+1)*250), 1] = as.numeric(d[1,]);
  th[ (ii*250+1):((ii+1)*250), 2] = as.numeric(d[2,]);
  th[ (ii*250+1):((ii+1)*250), 3] = as.numeric(d[3,]);
  ll[ (ii*250+1):((ii+1)*250), 1] = as.numeric(d[4,]);
}

llmax = max(ll)
w = exp(ll - llmax)
w = w / sum( w )

layout(matrix(1:3, 3, 1, byrow = FALSE))  
par(mar=c(4,5,1,1))

plot(density(th[,1],weights=w,to=1.0),xlab=expression(phi),main="",col="darkgreen",lwd=3,xlim=c(0.9,1),ylim=c(-0.05,15))
points(th[,1],rep(-0.05,length(th[,1])),pch=19,col="darkgrey",cex=0.2)
lines(c(-1,1),c(1,1),lwd=3,col="steelblue")
abline(v=sum(w * th[,1]),lty="dotted")

plot(density(th[,2],weights=w,from=0.0),xlab=expression(sigma[v]),main="",col="darkred",lwd=3,xlim=c(0,0.25),ylim=c(-0.05,15))
points(as.numeric(th[,2]),rep(-0.05,length(th[,2])),pch=19,col="darkgrey",cex=0.2)
xx = seq(0,0.5,0.01)
yy = dgamma(xx,shape=2,scale=0.1)
lines(xx,yy,lwd=3,col="steelblue")
abline(v=sum(w * th[,2]),lty="dotted")

plot(density(th[,3],weights=w,to=1.0,from=0.0),xlab=expression(beta),main="",col="tan3",lwd=3,xlim=c(0.55,0.95),ylim=c(-0.05,6))
points(as.numeric(th[,3]),rep(-0.05,length(th[,3])),pch=19,col="darkgrey",cex=0.2)
xx = seq(0,1,0.01)
yy = dgamma(xx,shape=7,scale=0.1)
lines(xx,yy,lwd=3,col="steelblue")
abline(v=sum(w * th[,3]),lty="dotted")

print("Parameter estimates")
sum(w * th[,1])
sum(w * th[,2])
sum(w * th[,3])
