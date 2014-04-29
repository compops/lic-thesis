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
# GPO using different acqusition rules
# Example 4.7 in Section 4.4.2 
#
########################################################################

import numpy as np
import scipy as sp
import pylab as pb
import GPy
import pandas

########################################################################
# Create the underlying function and evaluate it on a grid
########################################################################
x = np.linspace(0.,10.,500); x = x[:,None];
y = (3*np.cos(x/2) + np.sin(2+x/4) )**2;

# Select some random samples of the function and add some noise to the observations
xred = x[[ 20, 197, 317, 149,  89, 443,  53, 353, 395, 186, 292, 359, 326, 369,  99]]
yred = y[[ 20, 197, 317, 149,  89, 443,  53, 353, 395, 186, 292, 359, 326, 369,  99]] + np.sqrt(2) * np.random.normal(size=(15,1))

# Calculate the GPO model using all the data
k3 = GPy.kern.bias(1) + GPy.kern.Matern52(1, lengthscale=1);     
m3 = GPy.models.GPRegression(xred,yred,k3); m3.optimize();

########################################################################
# Fit the GP to the first 3 observations and predict the values of the grid
########################################################################

m3 = GPy.models.GPRegression(xred[0:3],yred[0:3],k3);
Mup, ys2, up95, lo95 = m3.predict( np.array(x) ); 

# Compute the values of the different AQs over the grid
mumax = max(Mup);
Z     = ( Mup - mumax - 0.01 ) / np.sqrt( ys2 )
PI    = sp.stats.norm.cdf( Z , 0, 1)
EI    = Z * sp.stats.norm.cdf( Z, 0, 1) + np.sqrt( ys2 ) * sp.stats.norm.pdf( Z, 0, 1)
UCB   = up95;

# Export data to R for plotting
out = np.hstack((Mup,up95,lo95,PI,EI,UCB));
pandas.DataFrame(out).to_csv("../figures/ch4-example-gp-acqfunc/ch4-example-gp-acqfunc-3d.csv");

########################################################################
# Fit the GP to the first 5 observations and predict the values of the grid
########################################################################

m3 = GPy.models.GPRegression(xred[0:5],yred[0:5],k3);
Mup, ys2, up95, lo95 = m3.predict( np.array(x) ); 

# Compute the values of the different AQs over the grid
mumax = max(Mup);
Z     = ( Mup - mumax - 0.01 ) / np.sqrt( ys2 )
PI    = sp.stats.norm.cdf( Z , 0, 1)
EI    = Z * sp.stats.norm.cdf( Z, 0, 1) + np.sqrt( ys2 ) * sp.stats.norm.pdf( Z, 0, 1)
UCB   = up95;

# Export data to R for plotting
out = np.hstack((Mup,up95,lo95,PI,EI,UCB));
pandas.DataFrame(out).to_csv("../figures/ch4-example-gp-acqfunc/ch4-example-gp-acqfunc-5d.csv");

########################################################################
# Fit the GP to the first 10 observations and predict the values of the grid
########################################################################

m3 = GPy.models.GPRegression(xred[0:10],yred[0:10],k3);
Mup, ys2, up95, lo95 = m3.predict( np.array(x) ); 

# Compute the values of the different AQs over the grid
mumax = max(Mup);
Z     = ( Mup - mumax - 0.01 ) / np.sqrt( ys2 )
PI    = sp.stats.norm.cdf( Z , 0, 1)
EI    = Z * sp.stats.norm.cdf( Z, 0, 1) + np.sqrt( ys2 ) * sp.stats.norm.pdf( Z, 0, 1)
UCB   = up95;

# Export data to R for plotting
out = np.hstack((Mup,up95,lo95,PI,EI,UCB));
pandas.DataFrame(out).to_csv("../figures/ch4-example-gp-acqfunc/ch4-example-gp-acqfunc-10d.csv");

########################################################################
# Export data to R for plotting
########################################################################

out = np.hstack((x,y));
pandas.DataFrame(out).to_csv("../figures/ch4-example-gp-acqfunc/ch4-example-gp-acqfunc-data.csv");

out = np.hstack((xred,yred));
pandas.DataFrame(out).to_csv("../figures/ch4-example-gp-acqfunc/ch4-example-gp-acqfunc-datared.csv");

########################################################################
# End of file
########################################################################