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
# Log-likelihood, score and information matrix in the LGSS model
# Example 2.6 in Section 2.3
#
########################################################################

d <- read.table("ch2-example-likelihoodtheory.csv",sep=",",header=F)
d=d[-1,];d=d[,-1];

layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))  
par(mar=c(4,5,1,1)) 

plot(d[,1],d[,2],type="l",xlab=expression(phi),ylab="Log-likelihood function",lwd=3,col="steelblue",ylim=c(-400,-340),cex.axis=1); 
abline(v=0.5,lty="dotted");

plot(d[,1],d[,3],type="l",xlab=expression(phi),ylab="Score function",lwd=3,col="darkred",ylim=c(-200,200),cex.axis=1); 
abline(v=0.5,lty="dotted"); abline(h=0,lty="dotted");

plot(d[,1],d[,4],type="l",xlab=expression(phi),ylab="Expected information",lwd=3,col="darkgreen",ylim=c(0,1000),yaxt="n",cex.axis=1); 
axis(2,at=seq(100,1000,100),cex.axis=1)
abline(v=0.5,lty="dotted"); abline(h=0,lty="dotted");

