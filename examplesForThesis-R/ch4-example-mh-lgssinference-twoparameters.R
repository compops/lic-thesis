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
# Parameter inference in the LGSS model
# Example 4.1 in Section 4.2
#
########################################################################

d <- read.table("ch4-example-lgssinference.csv",header=T,sep=",",stringsAsFactors=F)[,-1]

layout(matrix(1:4, 2, 2, byrow = TRUE))  
par(mar=c(4,5,1,1)) 

plot(as.numeric(d[1:500,1]),type="l",ylab=expression(phi),xlab="Iteration",cex.axis=1,lwd=2,col="darkgreen"); 
abline(h=0.50,lty="dotted");

hist(as.numeric(d[1000:10000,1]),breaks=floor(sqrt(9000)),main="",freq=F,xlab=expression(phi),border="darkgrey",col="darkgrey");
lines(density(as.numeric(d[1000:10000,1])),lwd=4,col="darkgreen")

plot(as.numeric(d[1:500,2]),type="l",ylab=expression(sigma[v]),xlab="Iteration",cex.axis=1,lwd=2,col="darkred"); 
abline(h=1.00,lty="dotted");

hist(as.numeric(d[1000:10000,2]),breaks=floor(sqrt(9000)),main="",freq=F,xlab=expression(sigma[v]),border="darkgrey",col="darkgrey");
lines(density(as.numeric(d[1000:10000,2])),lwd=4,col="darkred")

########################################################################
# End of file
########################################################################