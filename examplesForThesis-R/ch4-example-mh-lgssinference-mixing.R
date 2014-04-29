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
# Parameter inference in the LGSS model
# Example 4.1 in Section 4.2
#
########################################################################

d1 <- read.table("ch4-example-lgssinference-mixing-small.csv",header=T,sep=",")[,-1]
d2 <- read.table("ch4-example-lgssinference-mixing-right.csv",header=T,sep=",")[,-1]
d3 <- read.table("ch4-example-lgssinference-mixing-large.csv",header=T,sep=",")[,-1]

layout(matrix(1:6, 3, 2, byrow = TRUE))  
par(mar=c(4,5,1,1)) 

plot(d1[1:500,1],type="l",ylab=expression(phi),xlab="Iteration",cex.axis=1,lwd=2,col="darkgreen"); 
text(520,0.12,expression(epsilon * " =0.01"),col="darkgreen",pos=2)
abline(h=0.50,lty="dotted"); text(9500,0.12,expression(epsilon *"=0.01"))
acf(d1[1000:10000,1],lag.max=50,main="",ylab=expression("ACF of " * phi), ci.col="steelblue", lwd=3,col="tan3")

plot(d2[1:500,1],type="l",ylab=expression(phi),xlab="Iteration",cex.axis=1,lwd=2,col="darkred");
text(520,0.12,expression(epsilon * " =0.10"),col="darkred",pos=2)
abline(h=0.50,lty="dotted"); text(9500,0.12,expression(epsilon *"=0.10"))
acf(d2[1000:10000,1],lag.max=50,main="",ylab=expression("ACF of " * phi), ci.col="steelblue", lwd=3,col="tan3")

plot(d3[1:500,1],type="l",ylab=expression(phi),xlab="Iteration",cex.axis=1,lwd=2,col="steelblue");
text(520,0.12,expression(epsilon * " =1.00"),col="steelblue",pos=2)
abline(h=0.50,lty="dotted"); text(9500,0.12,expression(epsilon *"=1.00"))
acf(d3[1000:10000,1],lag.max=50,main="",ylab=expression("ACF of " * phi), ci.col="steelblue", lwd=3,col="tan3")

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
ess(d1[1000:10000,1],0,0)
ess(d2[1000:10000,1],0,0)
ess(d3[1000:10000,1],0,0)

########################################################################
# End of file
########################################################################