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
# GPO for ML inference in the GARCH(1,1) model
# Example 4.8 in Section 4.4.3
#
########################################################################

grid=seq(0.01,0.25,0.01)
toplott = c(3:7)

layout(matrix(1:15, 5, 3, byrow = FALSE))  
par(mar=c(4,5,1,3.5)) 

########################################################################
# PI rule
########################################################################

d <- read.table("ch4-example-gpo-garchinference-pi.csv",header=T,sep=",",stringsAsFactors=F)[,-1]
d2 <- read.table("ch4-example-gpo-garchinference-pi2.csv",header=T,sep=",",stringsAsFactors=F)[,-1]
d2 = as.numeric(d2)
d2th = d2[1:11];
d2ll = d2[12:21];

M = d[,1:10];
L = d[,11:20];
U = d[,21:30];
A = d[,30:40];

for ( ii in toplott ) {
  plot(grid,M[,ii],type="l",lwd=3,col="darkgreen",ylim=c(-300,-210),xlab=expression(alpha),ylab="Log-likelihood",xlim=c(0,0.25))
  polygon( c(grid,rev(grid)), c(U[,ii], rev(L[,ii])), col="lightgrey", border=NA)
  lines(grid,M[,ii],lwd=3,col="darkgreen")
  points(d2th[1:ii-1],d2ll[1:ii-1],pch=19,cex=1.5)
  
  par(new=T)
  AA = -as.numeric( A[,ii] )
  plot(grid,AA,lwd=2,col="darkgreen",lty="dotted",axes=F,xlab="",ylab="", ylim=c(0,2*max(AA)),type="l")
  axis(4, ylim=c(0,2*max(AA)),lwd=1)
  mtext(4,text="PI",line=2,cex=0.65)
}

########################################################################
# EI rule
########################################################################

d <- read.table("ch4-example-gpo-garchinference-ei.csv",header=T,sep=",",stringsAsFactors=F)[,-1]
d2 <- read.table("ch4-example-gpo-garchinference-ei2.csv",header=T,sep=",",stringsAsFactors=F)[,-1]
d2 = as.numeric(d2)
d2th = d2[1:11];
d2ll = d2[12:21];

M = d[,1:10];
L = d[,11:20];
U = d[,21:30];
A = d[,30:40];

for ( ii in toplott ) {
  plot(grid,M[,ii],type="l",lwd=3,col="darkred",ylim=c(-240,-210),xlab=expression(alpha),ylab="Log-likelihood",xlim=c(0,0.25))
  polygon( c(grid,rev(grid)), c(U[,ii], rev(L[,ii])), col="lightgrey", border=NA)
  lines(grid,M[,ii],lwd=3,col="darkred")
  points(d2th[1:ii-1],d2ll[1:ii-1],pch=19,cex=1.5)
  
  par(new=T)
  AA = -as.numeric( A[,ii] )
  plot(grid,AA,lwd=2,col="darkred",lty="dotted",axes=F,xlab="",ylab="", ylim=c(0,2*max(AA)),type="l")
  axis(4, ylim=c(0,2*max(AA)),lwd=1)
  mtext(4,text="EI",line=2,cex=0.65)
}

########################################################################
# UCB rule
########################################################################

d <- read.table("ch4-example-gpo-garchinference-ucb.csv",header=T,sep=",",stringsAsFactors=F)[,-1]
d2 <- read.table("ch4-example-gpo-garchinference-ucb2.csv",header=T,sep=",",stringsAsFactors=F)[,-1]
d2 = as.numeric(d2)
d2th = d2[1:11];
d2ll = d2[12:21];

M = d[,1:10];
L = d[,11:20];
U = d[,21:30];
A = d[,30:40];

for ( ii in toplott ) {
  plot(grid,M[,ii],type="l",lwd=3,col="steelblue",ylim=c(-300,-210),xlab=expression(alpha),ylab="Log-likelihood",xlim=c(0,0.25))
  polygon( c(grid,rev(grid)), c(U[,ii], rev(L[,ii])), col="lightgrey", border=NA)
  lines(grid,M[,ii],lwd=3,col="steelblue")
  points(d2th[1:ii-1],d2ll[1:ii-1],pch=19,cex=1.5)
  
  par(new=T)
  AA = -as.numeric( A[,ii] )
  plot(grid,AA,lwd=2,col="steelblue",lty="dotted",axes=F,xlab="",ylab="", ylim=c(1.3*max(AA),-150),type="l")
  axis(4, ylim=c(2*max(AA),0),lwd=1)
  mtext(4,text="UCB",line=2,cex=0.65)
}

########################################################################
# End of file
########################################################################