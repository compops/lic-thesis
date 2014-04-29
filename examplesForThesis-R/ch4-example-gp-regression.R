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
# GP regression
# Example 4.6 in Section 4.4.1
#
########################################################################

dat  = read.table("ch4-example-gp-regression-data.csv",sep=",",header=T)[,-1]
datr = read.table("ch4-example-gp-regression-datared.csv",sep=",",header=T)[,-1]
rbfa  = read.table("ch4-example-gp-regression-rbf-10d.csv",sep=",",header=T)[,-1]
m32a  = read.table("ch4-example-gp-regression-m32-10d.csv",sep=",",header=T)[,-1]
m52a  = read.table("ch4-example-gp-regression-m52-10d.csv",sep=",",header=T)[,-1]
rbfb  = read.table("ch4-example-gp-regression-rbf-20d.csv",sep=",",header=T)[,-1]
m32b  = read.table("ch4-example-gp-regression-m32-20d.csv",sep=",",header=T)[,-1]
m52b  = read.table("ch4-example-gp-regression-m52-20d.csv",sep=",",header=T)[,-1]

layout(matrix(1:6, 3, 2, byrow = FALSE)); par(mar=c(4,5,1,1)) 

plot(dat[,1],rbfa[,1],type="l",lwd=3,ylim=c(-8,30),xlim=c(0.30,9.65),xlab="x",ylab="f(x)",col="darkgreen"); 
polygon( c(dat[,1],rev(dat[,1])), c(rbfa[,3], rev(rbfa[,2])), col="lightgrey", border=NA)
lines(dat[,1],rbfa[,1],lwd=3,col="darkgreen"); lines(dat[,1],dat[,2],lwd=3,lty="dashed"); points(datr[1:5,1],datr[1:5,2],pch=19,cex=1.5)
text(10,30,"SE, N=5",pos=2);

plot(dat[,1],m52a[,1],type="l",lwd=3,ylim=c(-8,30),xlim=c(0.30,9.65),col="darkred",xlab="x",ylab="f(x)"); 
polygon( c(dat[,1],rev(dat[,1])), c(m52a[,3], rev(m52a[,2])), col="lightgrey", border=NA)
lines(dat[,1],m52a[,1],lwd=3,col="darkred"); lines(dat[,1],dat[,2],lwd=3,lty="dashed"); points(datr[1:5,1],datr[1:5,2],pch=19,cex=1.5)
text(10,30,"Matérn 5/2, N=5",pos=2);

plot(dat[,1],m32a[,1],type="l",lwd=3,ylim=c(-8,30),xlim=c(0.30,9.65),col="steelblue",xlab="x",ylab="f(x)"); 
polygon( c(dat[,1],rev(dat[,1])), c(m32a[,3], rev(m32a[,2])), col="lightgrey", border=NA)
lines(dat[,1],m32a[,1],lwd=3,col="steelblue"); lines(dat[,1],dat[,2],lwd=3,lty="dashed"); points(datr[1:5,1],datr[1:5,2],pch=19,cex=1.5)
text(10,30,"Matérn 3/2, N=5",pos=2);

plot(dat[,1],rbfb[,1],type="l",lwd=3,ylim=c(-8,30),xlim=c(0.30,9.65),xlab="x",ylab="f(x)",col="darkgreen"); 
polygon( c(dat[,1],rev(dat[,1])), c(rbfb[,3], rev(rbfb[,2])), col="lightgrey", border=NA)
lines(dat[,1],rbfb[,1],lwd=3,col="darkgreen"); lines(dat[,1],dat[,2],lwd=3,lty="dashed"); points(datr[1:15,1],datr[1:15,2],pch=19,cex=1.5)
text(10,30,"SE, N=15",pos=2);

plot(dat[,1],m52b[,1],type="l",lwd=3,ylim=c(-8,30),xlim=c(0.30,9.65),col="darkred",xlab="x",ylab="f(x)"); 
polygon( c(dat[,1],rev(dat[,1])), c(m52b[,3], rev(m52b[,2])), col="lightgrey", border=NA)
lines(dat[,1],m52b[,1],lwd=3,col="darkred"); lines(dat[,1],dat[,2],lwd=3,lty="dashed"); points(datr[1:15,1],datr[1:15,2],pch=19,cex=1.5)
text(10,30,"Matérn 5/2, N=15",pos=2);

plot(dat[,1],m32b[,1],type="l",lwd=3,ylim=c(-8,30),xlim=c(0.30,9.65),col="steelblue",xlab="x",ylab="f(x)"); 
polygon( c(dat[,1],rev(dat[,1])), c(m32b[,3], rev(m32b[,2])), col="lightgrey", border=NA)
lines(dat[,1],m32b[,1],lwd=3,col="steelblue"); lines(dat[,1],dat[,2],lwd=3,lty="dashed"); points(datr[1:15,1],datr[1:15,2],pch=19,cex=1.5)
text(10,30,"Matérn 3/2, N=15",pos=2);

########################################################################
# End of file
########################################################################