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
# IS for Bayesian parameter inference in the LGSS model
# Example 3.1 in Section 3.2 
#
########################################################################

d <- read.table("ch3-example-importancesampling-lgss.csv",header=T,sep=",",stringsAsFactors=F)[,-1]

layout(matrix(1:2, 1, 2, byrow = FALSE))  
par(mar=c(4,5,1,1))

plot(density(as.numeric(d[1,]),weights=as.numeric(d[3,])),xlab=expression(phi),main="",col="darkgreen",lwd=3,xlim=c(-1,1),ylim=c(-0.05,3))
points(as.numeric(d[1,]),rep(-0.05,length(d[1,])),pch=19,col="darkgrey",cex=0.2)
lines(c(-1,1),c(1,1),lwd=3,col="steelblue")
abline(v=0.5,lty="dotted")
plot(density(as.numeric(d[2,]),weights=as.numeric(d[3,])),xlab=expression(sigma[v]),main="",col="darkred",lwd=3,xlim=c(0,4),ylim=c(-0.05,2))
points(as.numeric(d[2,]),rep(-0.05,length(d[1,])),pch=19,col="darkgrey",cex=0.2)
abline(v=1.0,lty="dotted")

xx = seq(0,5,0.01)
yy = dgamma(xx,shape=1,scale=1)
lines(xx,yy,lwd=3,col="steelblue")

########################################################################
# End of file
########################################################################