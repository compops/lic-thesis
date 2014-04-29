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
# Input design in the LGSS model using GPO
# Example 4.10 in Section 4.4.3
#
########################################################################

d <- read.table("ch4-example-lgss-inputdesign.csv",sep=",",header=T,stringsAsFactors=F)[,-1]
th <- read.table("ch4-example-lgss-inputdesign-th.csv",sep=",",header=T,stringsAsFactors=F)[,-1]
d = t(as.numeric(d))

# The parameter estimates (may differ between runs)
alpha1th = 0.23
alpha2th = 0.19

layout(matrix(1:3, 3, 1, byrow = TRUE))  
par(mar=c(4,5,1,2))  

xs = seq(0, 1.00, length.out=101)
ys = seq(0, 1.00, length.out=101)
contour(ys,xs,matrix(-d, nrow=101, ncol=101, byrow=T),nlevels = 30, lwd=2, 
        xlab = expression(alpha[1]), ylab = expression(alpha[2]), col=rainbow(60,start=0,end=0.9),
        drawlabels = F, xlim=c(0.08,0.92),ylim=c(0.08,0.92),xaxt="n",yaxt="n")
axis(1,at=seq(0,1,0.20))
axis(2,at=seq(0,1,0.20))

points(th[,1],th[,2],pch=19,cex=0.8)
abline(v=0.09,lty="dotted",lwd=2)
abline(h=0.30,lty="dotted",lwd=2)

plot(th[,1],ylab="Sample point",xlab="Iteration",type="l",ylim=c(0,1),lwd=2)
text(50,0.4,expression(alpha[1]),pos=1)
lines(th[,2],col="darkred",lwd=2); 
text(50,0.15,expression(alpha[2]),pos=1,col="darkred")

# Generate a realisation of the input
u1 = runif(100,min=0,max=1)
u2 = runif(100,min=-alpha2th,max=1-alpha2th)
u = (alpha1th < u1) * sign(u2)

plot(u,type="s",lwd=2,col="darkgreen",xlab="Time",ylab="Input")

########################################################################
# End of file
########################################################################
