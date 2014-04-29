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
# GPO using different acquisition rules
# Example 4.7 in Section 4.4.2 
#
########################################################################

# Helper function to normalise the AQs
normalise <- function(x){
  return( x/max(x) - min(x)/max(x) );
}


dat   = read.table("ch4-example-gp-acqfunc-data.csv",sep=",",header=T)[,-1]
datr  = read.table("ch4-example-gp-acqfunc-datared.csv",sep=",",header=T)[,-1]
gp3   = read.table("ch4-example-gp-acqfunc-3d.csv",sep=",",header=T)[,-1]
gp5   = read.table("ch4-example-gp-acqfunc-5d.csv",sep=",",header=T)[,-1]
gp10  = read.table("ch4-example-gp-acqfunc-10d.csv",sep=",",header=T)[,-1]

layout(matrix(1:6, 3, 2, byrow = FALSE)); par(mar=c(4,5,1,1)) 

plot(dat[,1],gp3[,1],type="l",lwd=3,ylim=c(-8,30),xlim=c(0.30,9.65),xlab="x",ylab="f(x)"); 
polygon( c(dat[,1],rev(dat[,1])), c(gp3[,3], rev(gp3[,2])), col="lightgrey", border=NA)
lines(dat[,1],gp3[,1],lwd=3); lines(dat[,1],dat[,2],lwd=3,lty="dashed"); points(datr[1:3,1],datr[1:3,2],pch=19,cex=1.5)
text(10,30,"Matérn 5/2, N=3",pos=2);

plot(dat[,1],gp5[,1],type="l",lwd=3,ylim=c(-8,30),xlim=c(0.30,9.65),col="black",xlab="x",ylab="f(x)");
polygon( c(dat[,1],rev(dat[,1])), c(gp5[,3], rev(gp5[,2])), col="lightgrey", border=NA)
lines(dat[,1],gp5[,1],lwd=3,col="black"); lines(dat[,1],dat[,2],lwd=3,lty="dashed"); points(datr[1:5,1],datr[1:5,2],pch=19,cex=1.5)
text(10,30,"Matérn 5/2, N=5",pos=2);

plot(dat[,1],gp10[,1],type="l",lwd=3,ylim=c(-8,30),xlim=c(0.30,9.65),col="black",xlab="x",ylab="f(x)");
polygon( c(dat[,1],rev(dat[,1])), c(gp10[,3], rev(gp10[,2])), col="lightgrey", border=NA)
lines(dat[,1],gp10[,1],lwd=3,col="black"); lines(dat[,1],dat[,2],lwd=3,lty="dashed"); points(datr[1:10,1],datr[1:10,2],pch=19,cex=1.5)
text(10,30,"Matérn 5/2, N=10",pos=2);

plot(dat[,1],normalise(gp3[,3]),type="l",lwd=3,ylim=c(0,1.2),xlim=c(0.30,9.65),xlab="x",ylab="AQ(x)",col="darkgreen"); 
lines(dat[,1],normalise(gp3[,4]),lwd=3,col="darkred")
lines(dat[,1],normalise(gp3[,5]),lwd=3,col="steelblue")
text(10,0.85,"PI",pos=2,font=2,cex=1.2,col="darkgreen");
text(10,1.00,"UCB",pos=2,col="steelblue",font=2,cex=1.2);
text(10,0.35,"EI",pos=2,col="darkred",font=2,cex=1.2);

plot(dat[,1],normalise(gp5[,3]),type="l",lwd=3,ylim=c(0,1.2),xlim=c(0.30,9.65),xlab="x",ylab="AQ(x)",col="darkgreen");
lines(dat[,1],normalise(gp5[,4])-0.02,lwd=3,col="darkred")
lines(dat[,1],normalise(gp5[,5]),lwd=3,col="steelblue")

plot(dat[,1],normalise(gp10[,3]),type="l",lwd=3,ylim=c(0,1.2),xlim=c(0.30,9.65),xlab="x",ylab="AQ(x)",col="darkgreen");
lines(dat[,1],normalise(gp10[,4])-0.03,lwd=3,col="darkred")
lines(dat[,1],normalise(gp10[,5]),lwd=3,col="steelblue")

########################################################################
# End of file
########################################################################