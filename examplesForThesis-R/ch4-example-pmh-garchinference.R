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
# PMH0 for parameter inference in the GARCH(1,1) model
# Example 4.4 in Section 4.3
#
########################################################################

d <- read.table("pmh0-run.csv",header=T,sep=",",stringsAsFactors=F)[,-1]

layout(matrix(c(1,1,2,3,3,4,5,5,6,7,7,8), 4, 3, byrow = TRUE))  
par(mar=c(4,5,1,1)) 

plot(seq(5000,10000),d[5000:10000,1],col="darkgreen",type="l",ylab=expression(alpha),xlab="Iteration",cex.axis=1,lwd=2); 
abline(h=mean(d[5000:20000,1]),lty="dotted");

hist(d[5000:20000,1],breaks=floor(sqrt(15000)),main="",freq=F,xlab=expression(phi),border="darkgrey",col="darkgrey",xlim=c(0,0.020));
lines(density(d[5000:20000,1]),lwd=4,col="darkgreen")

plot(seq(5000,10000),d[5000:10000,2],,col="steelblue",type="l",ylab=expression(beta),xlab="Iteration",cex.axis=1,lwd=2); 
abline(h=mean(d[5000:20000,2]),lty="dotted");

hist(d[5000:20000,2],breaks=floor(sqrt(15000)),main="",freq=F,xlab=expression(beta),border="darkgrey",col="darkgrey",xlim=c(0.10,0.20));
lines(density(d[5000:20000,2]),lwd=4,col="steelblue")

plot(seq(5000,10000),d[5000:10000,3],,col="darkred",type="l",ylab=expression(gamma),xlab="Iteration",cex.axis=1,lwd=2); 
abline(h=mean(d[5000:20000,3]),lty="dotted");

hist(d[5000:20000,3],breaks=floor(sqrt(15000)),main="",freq=F,xlab=expression(gamma),border="darkgrey",col="darkgrey",xlim=c(0.80,0.90));
lines(density(d[5000:20000,3]),lwd=4,col="darkred")

plot(seq(5000,10000),d[5000:10000,4],col="tan3",type="l",ylab=expression(tau),xlab="Iteration",cex.axis=1,lwd=2); 
abline(h=mean(d[5000:20000,4]),lty="dotted");

hist(d[5000:20000,4],breaks=floor(sqrt(15000)),main="",freq=F,xlab=expression(tau),border="darkgrey",col="darkgrey",xlim=c(0.55,0.70));
lines(density(d[5000:20000,4]),lwd=4,col="tan3")

print("Posterior mean for each parameter:")
mean(d$th0[5000:20000])
mean(d$th1[5000:20000])
mean(d$th2[5000:20000])
mean(d$th3[5000:20000])

print("Acceptance rate:")
mean(d$X14[5000:20000])

########################################################################
# Calculate the effective sample size
########################################################################

ess <- function(x,xt,lag=0) {
  x = x - xt;
  xcorr = acf(x, lag.max = floor( length(x)/2 ), plot=F )
  lim   = 2/sqrt(length(x));
  
  if ( sum(is.nan(xcorr$acf)) > 0 ) {
    out = NA;
  } else {
    if (lag==0) {
      newli = min(which(xcorr$acf < lim)[1],floor( length(x) / 10),na.rm=T);
    } else {
      newli = lag
    }
    out = 1 + 2 * sum( xcorr$acf[1:newli] );
  }
  S = length(x);
  return(S/out);
}

print("ESS for each parameter:")
ess(d[5000:20000,1],0,0)
ess(d[5000:20000,2],0,0)
ess(d[5000:20000,3],0,0)
ess(d[5000:20000,4],0,0)

########################################################################
# End of file
########################################################################


