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
# State inference in the earthquake count model
# Example 3.2 in Section 3.3.2
#
########################################################################

d <- read.table("ch3-example-earthquakefiltering.csv",header=F,sep=",",stringsAsFactors=FALSE)[,-1]
d2 <- read.table("earthquakes.txt",sep="\t")
d3 <- read.table("ch3-example-earthquakefiltering-pred.csv",header=F,sep=",",stringsAsFactors=FALSE)[,-1]
d <- d[-1,]
d3 <- d3[-1,]

########################################################################
# Compute predictions
########################################################################

nPred = 10;
mpred = matrix(0,nrow=nPred,ncol=1)
spred = matrix(0,nrow=nPred,ncol=1)
meart = matrix(0,nrow=nPred,ncol=1)
peart = matrix(0,nrow=nPred,ncol=1)
for (tt in 1:nPred) {
  mpred[tt] = mean(as.numeric(d3[,tt]));
  spred[tt] = sd(as.numeric(d3[,tt]));
  meart[tt] = mean( 17.64746228 * exp(as.numeric(d3[,tt])  ))
  peart[tt] = sd( 17.64746228 * exp(as.numeric(d3[,tt])  ))
}

########################################################################
# Plotting
########################################################################

layout(matrix(1:2, 2, 1, byrow = TRUE))  
par(mar=c(4,5,1,1)) 

plot(d2[,1],d2[,2],col="black",pch=19,cex=0.5,xlab="Year",ylab="Number of major earthquakes",cex.lab=0.8,cex.axis=0.8,xlim=c(1900,2013+nPred),ylim=c(0,50))
lines(d2[,1],as.numeric(d[,2]),lwd=3,col="darkred")
polygon( c(seq(2013,(2013+nPred-1),1),rev(seq(2013,(2013+nPred-1),1))), c( meart-1.96*peart, rev(meart+1.96*peart)), col="lightgrey", border=NA)
lines(seq(2013,(2013+nPred-1),1),meart,lwd=3,col="darkgreen")

plot(d2[,1],d[,1],type="l",col="steelblue",xlab="Year",ylab="Estimated latent state",lwd=3,cex.lab=0.8,cex.axis=0.8,xlim=c(1900,2013+nPred),ylim=c(-1,1))

polygon( c(seq(2013,(2013+nPred-1),1),rev(seq(2013,(2013+nPred-1),1))), c(mpred-1.96*spred, rev(mpred+1.96*spred)), col="lightgrey", border=NA)
lines(seq(2013,(2013+nPred-1),1),mpred,lwd=3,col="darkgreen")

########################################################################
# End of file
########################################################################
