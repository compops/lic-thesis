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
# Bias and variance of the log-likelihood estimate
# Example 3.4 in Section 3.3.4
#
########################################################################

d <- read.table("ch3-example-llestimationCLT.csv",sep=",",header=F)
d=d[-1,];d=d[,-1];

layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))  
par(mar=c(4,5,1,1)) 

hist(d,breaks=floor(sqrt(1000)),main="",freq=F,xlab="Error in the log-likelihood estimate",xlim=c(-1,1),col="grey",ylim=c(0,2),border="grey")

x   <- seq(-1,1,length=1000)
y   <- dnorm(x,mean(d), sd(d))
lines(x,y,type="l",lwd=5,col="steelblue")

boxplot(d,ylab="Error in the log-likelihood estimate",horizontal=F)
qqnorm(d,main=""); qqline(d,lwd=3,col="steelblue");

var(d)
mean(d)

########################################################################
# End of file
########################################################################