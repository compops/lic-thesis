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
# Mixing property in the LGSS model
# Example 3.5 in Section 3.4.1
#
########################################################################

stateseqs1  <- read.table("ch3-example-forgettingproperties1.csv",header=T,sep=",",stringsAsFactors=F)[,-1]
stateseqs2  <- read.table("ch3-example-forgettingproperties2.csv",header=T,sep=",",stringsAsFactors=F)[,-1]
stateseqs3  <- read.table("ch3-example-forgettingproperties3.csv",header=T,sep=",",stringsAsFactors=F)[,-1]

T = 21;
N = 10;

layout(matrix(1:3, 3, 1, byrow = TRUE))  
par(mar=c(4,5,1,1)) 

plot(0:14,stateseqs1[,1],xlim=c(0,8),ylim=c(-60,60),xlab="Time",ylab="State",cex.axis=0.95,type="l",lwd=2,col="steelblue")

for (ii in 2:20) {
  lines(0:14,stateseqs1[,ii],lwd=2,col="steelblue")
}
text(8.2,55,expression(phi * "=0.2"),col="steelblue",font=2,pos=2,cex=1.4)

plot(0:14,stateseqs2[,1],xlim=c(0,8),ylim=c(-60,60),xlab="Time",ylab="State",cex.axis=0.95,type="l",lwd=2,col="darkred")

for (ii in 2:20) {
  lines(0:14,stateseqs2[,ii],lwd=2,col="darkred")
}
text(8.2,55,expression(phi * "=0.5"),col="darkred",font=2,pos=2,cex=1.4)

plot(0:14,stateseqs3[,1],xlim=c(0,8),ylim=c(-60,60),xlab="Time",ylab="State",cex.axis=0.95,type="l",lwd=2,col="darkgreen")

for (ii in 2:20) {
  lines(0:14,stateseqs3[,ii],lwd=2,col="darkgreen")
}
text(8.2,55,expression(phi * "=0.8"),col="darkgreen",font=2,pos=2,cex=1.4)

########################################################################
# End of file
########################################################################
