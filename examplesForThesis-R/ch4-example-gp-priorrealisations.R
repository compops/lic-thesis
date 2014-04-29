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
# GP kernels
# Example 4.5 in Section 4.4.1
#
########################################################################

dx = read.table("ch4-example-gp-priorrealisations-x.csv",sep=",",header=T)[,-1]
dz1 = read.table("ch4-example-gp-priorrealisations-z1.csv",sep=",",header=T)[,-1]
dz2 = read.table("ch4-example-gp-priorrealisations-z2.csv",sep=",",header=T)[,-1]
dz3 = read.table("ch4-example-gp-priorrealisations-z3.csv",sep=",",header=T)[,-1]
dz4 = read.table("ch4-example-gp-priorrealisations-z4.csv",sep=",",header=T)[,-1]
dz5 = read.table("ch4-example-gp-priorrealisations-z5.csv",sep=",",header=T)[,-1]
dz6 = read.table("ch4-example-gp-priorrealisations-z6.csv",sep=",",header=T)[,-1]

layout(matrix(1:6, 3, 2, byrow = FALSE)); par(mar=c(4,5,1,1))

plot(dx,dz1[,1],type="l",lwd=3,xlab="x",ylab="f(x)",ylim=c(-3,3),col="darkgreen"); lines(dx,dz1[,2],lwd=3,col="darkgreen"); lines(dx,dz1[,3],lwd=3,col="darkgreen"); text(10,3,"SE, l=1",pos=2);
plot(dx,dz3[,1],type="l",lwd=3,xlab="x",ylab="f(x)",ylim=c(-3,3),col="darkred"); lines(dx,dz3[,2],lwd=3,col="darkred"); lines(dx,dz3[,3],lwd=3,col="darkred"); text(10,3,"Matérn 5/2, l=1",pos=2);
plot(dx,dz2[,1],type="l",lwd=3,xlab="x",ylab="f(x)",ylim=c(-3,3),col="steelblue"); lines(dx,dz2[,2],lwd=3,col="steelblue"); lines(dx,dz2[,3],lwd=3,col="steelblue"); text(10,3,"Matérn 3/2, l=1",pos=2);
plot(dx,dz4[,1],type="l",lwd=3,xlab="x",ylab="f(x)",ylim=c(-3,3),col="darkgreen"); lines(dx,dz4[,2],lwd=3,col="darkgreen"); lines(dx,dz4[,3],lwd=3,col="darkgreen"); text(10,3,"SE, l=3",pos=2);
plot(dx,dz6[,1],type="l",lwd=3,xlab="x",ylab="f(x)",ylim=c(-3,3),col="darkred"); lines(dx,dz6[,5],lwd=3,col="darkred"); lines(dx,dz6[,6],lwd=3,col="darkred"); text(10,3,"Matérn 5/2, l=3",pos=2);
plot(dx,dz5[,1],type="l",lwd=3,xlab="x",ylab="f(x)",ylim=c(-3,3),col="steelblue"); lines(dx,dz5[,2],lwd=3,col="steelblue"); lines(dx,dz5[,3],lwd=3,col="steelblue"); text(10,3,"Matérn 3/2, l=3",pos=2);

########################################################################
# End of file
########################################################################