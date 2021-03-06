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

import numpy as np
import pylab as pb
import GPy
import pandas

########################################################################
# Create the underlying function and evaluate it on a grid
########################################################################
x = np.linspace(0.,10.,500); x = x[:,None];
y = (3*cos(x/2) + sin(2+x/4) )**2;

# Select some random samples of the function and add some noise to the observations
xred = x[[ 20, 197, 317, 149,  89, 443,  53, 353, 395, 186, 292, 359, 326, 369,  99]]
yred = y[[ 20, 197, 317, 149,  89, 443,  53, 353, 395, 186, 292, 359, 326, 369,  99]] + np.sqrt(2) * np.random.normal(size=(15,1))
plot(x,y); plot(xred,yred,'ro')

########################################################################
# Fit the GP to the first 5 observations using the SE kernel
########################################################################
k1 = GPy.kern.bias(d) + GPy.kern.rbf(d, lengthscale=1);      
m1 = GPy.models.GPRegression(xred[0:5],yred[0:5],k1); m1.optimize();

# Evaluate the predictive distribution on a grid
Mup, ys2, up95, lo95 = m1.predict( np.array(x) ); 

# Export data to R for plotting
out = hstack((Mup,up95,lo95));
pandas.DataFrame(out).to_csv("../figures/ch4-example-gp-regression/ch4-example-gp-regression-rbf-10d.csv");

########################################################################
# Fit the GP to the first 5 observations using the Matern 3/2 kernel
########################################################################
k2 = GPy.kern.bias(d) + GPy.kern.Matern32(d, lengthscale=1);     
m2 = GPy.models.GPRegression(xred[0:5],yred[0:5],k2); m2.optimize();

# Evaluate the predictive distribution on a grid
Mup, ys2, up95, lo95 = m2.predict( np.array(x) ); 

# Export data to R for plotting
out = hstack((Mup,up95,lo95));
pandas.DataFrame(out).to_csv("../figures/ch4-example-gp-regression/ch4-example-gp-regression-m32-10d.csv");

########################################################################
# Fit the GP to the first 5 observations using the Matern 5/2 kernel
########################################################################
k3 = GPy.kern.bias(d) + GPy.kern.Matern52(d, lengthscale=1);     
m3 = GPy.models.GPRegression(xred[0:5],yred[0:5],k3); m3.optimize();

# Evaluate the predictive distribution on a grid
Mup, ys2, up95, lo95 = m3.predict( np.array(x) ); 

# Export data to R for plotting
out = hstack((Mup,up95,lo95));
pandas.DataFrame(out).to_csv("../figures/ch4-example-gp-regression/ch4-example-gp-regression-m52-10d.csv");

########################################################################
# Fit the GP to the 15 observations using the SE kernel
########################################################################
k1 = GPy.kern.bias(d) + GPy.kern.rbf(d, lengthscale=1);      
m1 = GPy.models.GPRegression(xred,yred,k1); m1.optimize();

# Evaluate the predictive distribution on a grid
Mup, ys2, up95, lo95 = m1.predict( np.array(x) ); 

# Export data to R for plotting
out = hstack((Mup,up95,lo95));
pandas.DataFrame(out).to_csv("../figures/ch4-example-gp-regression/ch4-example-gp-regression-rbf-20d.csv");

########################################################################
# Fit the GP to the 15 observations using the Matern 3/2 kernel
########################################################################
k2 = GPy.kern.bias(d) + GPy.kern.Matern32(d, lengthscale=1);     
m2 = GPy.models.GPRegression(xred,yred,k2); m2.optimize();

# Evaluate the predictive distribution on a grid
Mup, ys2, up95, lo95 = m2.predict( np.array(x) ); 

# Export data to R for plotting
out = hstack((Mup,up95,lo95));
pandas.DataFrame(out).to_csv("../figures/ch4-example-gp-regression/ch4-example-gp-regression-m32-20d.csv");

########################################################################
# Fit the GP to the 15 observations using the Matern 5/2 kernel
########################################################################
k3 = GPy.kern.bias(d) + GPy.kern.Matern52(d, lengthscale=1);     
m3 = GPy.models.GPRegression(xred,yred,k3); m3.optimize();

# Evaluate the predictive distribution on a grid
Mup, ys2, up95, lo95 = m3.predict( np.array(x) ); 

# Export data to R for plotting
out = hstack((Mup,up95,lo95));
pandas.DataFrame(out).to_csv("../figures/ch4-example-gp-regression/ch4-example-gp-regression-m52-20d.csv");

# Export data to R for plotting
out = hstack((x,y));
pandas.DataFrame(out).to_csv("../figures/ch4-example-gp-regression/ch4-example-gp-regression-data.csv");

# Export data to R for plotting
out = hstack((xred,yred));
pandas.DataFrame(out).to_csv("../figures/ch4-example-gp-regression/ch4-example-gp-regression-datared.csv");

########################################################################
# End of file
########################################################################
