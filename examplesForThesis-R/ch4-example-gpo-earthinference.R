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
# GPO for ML inference in the earthquake count model
# Example 4.9 in Section 4.4.3 
#
########################################################################

d <- read.table("ch4-example-gpo-earthinference.csv",sep=",",header=T,stringsAsFactors=F)[,-1]

layout(matrix(1:4, 2, 2, byrow = TRUE))  
par(mar=c(4,5,1,1)) 
plot(50:200,as.numeric(d[50:200,1]),type="l",lwd=3,col="darkgreen",ylab=expression(hat(phi)),ylim=c(0,1),xlab="Iteration")
plot(50:200,as.numeric(d[50:200,2]),type="l",lwd=3,col="darkred",ylab=expression(hat(sigma)[v]),ylim=c(0,1),xlab="Iteration")
plot(50:200,as.numeric(d[50:200,3]),type="l",lwd=3,col="steelblue",ylab=expression(hat(beta)),ylim=c(10,20),xlab="Iteration")
plot(50:200,-as.numeric(d[50:200,4]),type="l",lwd=3,col="tan3",ylab="log-likelihood estimate",ylim=c(-420,-320),xlab="Iteration")

########################################################################
# End of file
########################################################################