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
# Score and information matrix in the Hull-White SV model
# Example 3.7 in Section 3.4.2
#
########################################################################

d <- read.table("ch3-example-hwsvllscoreinfo.csv",sep=",",header=F)
d=d[-1,];d=d[,-1];

layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))  
par(mar=c(4,5,1,1)) 

plot(d[,1],d[,2],type="l",xlab=expression(phi),ylab="Estimated log-likelihood function",lwd=3,col="steelblue",cex.axis=1); 
abline(v=0.95,lty="dotted");

plot(d[,1],d[,3],type="l",xlab=expression(phi),ylab="Estimated score function",lwd=3,col="darkred",ylim=c(-100,100),cex.axis=1); 
abline(v=0.95,lty="dotted"); abline(h=0,lty="dotted");

plot(d[,1],d[,4],type="l",xlab=expression(phi),ylab="Estimated observed information",lwd=3,col="darkgreen",cex.axis=1); 
abline(v=0.95,lty="dotted"); abline(h=0,lty="dotted");

########################################################################
# End of file
########################################################################