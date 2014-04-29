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
# Path degeneracy in the GARCH(1,1) model
# Example 3.3 in Section 3.3.2
#
########################################################################

bpfp  <- read.table("ch3-example-garchpathdegenercy-bpf10-p.csv",header=T,sep=",",stringsAsFactors=F)[,-1]
bpfa  <- read.table("ch3-example-garchpathdegenercy-bpf10-a.csv",header=T,sep=",",stringsAsFactors=F)[,-1]
bpfp2  <- read.table("ch3-example-garchpathdegenercy-bpf20-p.csv",header=T,sep=",",stringsAsFactors=F)[,-1]
bpfa2  <- read.table("ch3-example-garchpathdegenercy-bpf20-a.csv",header=T,sep=",",stringsAsFactors=F)[,-1]
fapfp <- read.table("ch3-example-garchpathdegenercy-fapf-p.csv",header=T,sep=",",stringsAsFactors=F)[,-1]
fapfa <- read.table("ch3-example-garchpathdegenercy-fapf-a.csv",header=T,sep=",",stringsAsFactors=F)[,-1]

# Number of time steps and particles to plot
T = 21;
N = 10;

layout(matrix(1:3, 3, 1, byrow = TRUE))  
par(mar=c(4,5,1,1)) 

########################################################################
# bPF with N = 10
########################################################################

plot(rep(T-1,N),bpfp[,T],xlim=c(0,20),ylim=c(-1.6,1.6),xlab="Time",ylab="State",cex.axis=1,pch=19,col="steelblue",cex=0.9)
for (ii in 1:N)
{
  att = ii;
  for (tt in seq(T-1,1,-1) )
  {
    at = bpfa[att,tt+1]+1;
    lines(c(tt,tt-1),c(bpfp[att,tt+1],bpfp[at,tt]),col="steelblue",lwd=2)
    att = at;
  }
}

for (ii in 1:N)
{
  att = ii;
  for (tt in seq(T-1,1,-1) )
  {
    at = bpfa[att,tt+1]+1;
    points(tt-1,bpfp[ii,tt],pch=19,col="grey35",cex=0.5)
    points(tt-1,bpfp[at,tt],pch=19,col="steelblue",cex=0.8)
    att = at;
  }
}

points(rep(T-1,N),bpfp[,T],pch=19,col="steelblue",cex=0.8)
text(-0.1,-1.5,"bPF, N=10",col="steelblue",font=2,pos=4)
abline(v=seq(0,20,5),lty="dotted")

########################################################################
# bPF with N = 20
########################################################################

N = 20;
plot(rep(T-1,N),bpfp2[,T],xlim=c(0,20),ylim=c(-1.6,1.6),xlab="Time",ylab="State",cex.axis=1,pch=19,col="darkred",cex=0.9)
for (ii in 1:N)
{
  att = ii;
  for (tt in seq(T-1,1,-1) )
  {
    at = bpfa2[att,tt+1]+1;
    lines(c(tt,tt-1),c(bpfp2[att,tt+1],bpfp2[at,tt]),col="darkred",lwd=2)
    att = at;
  }
}

for (ii in 1:N)
{
  att = ii;
  for (tt in seq(T-1,1,-1) )
  {
    at = bpfa2[att,tt+1]+1;
    points(tt-1,bpfp2[ii,tt],pch=19,col="grey35",cex=0.5)
    points(tt-1,bpfp2[at,tt],pch=19,col="darkred",cex=0.8)
    att = at;
  }
}

points(rep(T-1,N),bpfp2[,T],pch=19,col="darkred",cex=0.8)
text(-0.1,-1.5,"bPF, N=20",col="darkred",font=2,pos=4)
abline(v=seq(0,20,5),lty="dotted")

########################################################################
# faPF with N = 10
########################################################################

N = 10;
plot(rep(T-1,N),fapfp[,T],xlim=c(0,20),ylim=c(-1.6,1.6),xlab="Time",ylab="State",cex.axis=1,pch=19,col="darkgreen",cex=0.7)
for (ii in 1:N)
{
  att = ii;
  for (tt in seq(T-1,1,-1) )
  {
    at = fapfa[att,tt+1]+1;
    lines(c(tt,tt-1),c(fapfp[att,tt+1],fapfp[at,tt]),col="darkgreen",lwd=2)
    att = at;
  }
}

for (ii in 1:N)
{
  att = ii;
  for (tt in seq(T-1,1,-1) )
  {
    points(tt-1,fapfp[ii,tt],pch=19,col="grey35",cex=0.5)
    at = fapfa[att,tt+1]+1;
    points(tt-1,fapfp[at,tt],pch=19,col="darkgreen",cex=0.8)
    att = at;
  }
}
points(rep(T-1,N),fapfp[,T],pch=19,col="darkgreen",cex=0.8)
text(-0.1,-1.5,"faPF, N=10",col="darkgreen",font=2,pos=4)
abline(v=seq(0,20,5),lty="dotted")

########################################################################
# End of file
########################################################################