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
# State inference in the Hull-White SV model
# Example 3.6 in Section 3.4.1
#
########################################################################

d <- read.table("ch3-example-hwsvfilteringsmoothing.csv",header=F,sep=",",stringsAsFactors=F)[,-1]
d2 <- read.table("seOMXdata.csv",sep=",",stringsAsFactors=F)
d <- d[-1,]

layout(matrix(1:3, 3, 1, byrow = TRUE))  
par(mar=c(4,5,1,1)) 

plot(as.Date(d2[,1]),as.numeric(d[,5])/100,pch=19,cex=0.7,col="grey",ylim=c(-0.12,0.12),xlab="Year",ylab="Daily log-returns",xaxt="n")
lines(as.Date(d2[,1]),as.numeric(d[,1])/100,col="darkred",lwd=3)
lines(as.Date(d2[1:3565,1]),-as.numeric(d[1:3565,2])/100,col="steelblue",lwd=3)
r <- as.POSIXct(range(as.Date(d2[,1])), "years"); axis.Date(1, at = seq(r[1], r[2], by = "years"))

text(as.Date("2014-06-01"),0.1,"bPF",col="darkred",font=2,pos=2)
text(as.Date("2014-06-01"),-0.1,"FL smoother",col="steelblue",font=2,pos=2)

plot(as.Date(d2[,1]),as.numeric(d[,3])/100,type="l",col="darkred",ylim=c(-0.02,0.04),xlab="Year",ylab="Filter estimate of volatility",lwd=2,xaxt="n")
r <- as.POSIXct(range(as.Date(d2[,1])), "years"); axis.Date(1, at = seq(r[1], r[2], by = "years"))

plot(as.Date(d2[1:3565,1]),as.numeric(d[1:3565,4])/100,type="l",col="steelblue",ylim=c(-0.02,0.04),xlab="Year",ylab="Smoother estimate of volatility",lwd=2,xaxt="n")
r <- as.POSIXct(range(as.Date(d2[,1])), "years"); axis.Date(1, at = seq(r[1], r[2], by = "years"))

########################################################################
# End of file
########################################################################