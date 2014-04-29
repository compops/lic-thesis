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

import numpy as np
import pylab as pb
import GPy
import pandas

# Generate a grid to evaluate the kernels on
d  = 1;
X  = np.linspace(0.,10.,500); X = X[:,None];
mu = np.zeros((500)); 

# Use length scale 1 and generate a realisation
k = GPy.kern.rbf(d, lengthscale=1);      C = k.K(X,X); Z1 = np.random.multivariate_normal(mu,C,10)
k = GPy.kern.Matern32(d, lengthscale=1); C = k.K(X,X); Z2 = np.random.multivariate_normal(mu,C,10)
k = GPy.kern.Matern52(d, lengthscale=1); C = k.K(X,X); Z3 = np.random.multivariate_normal(mu,C,10)

# Use length scale 2 and generate a realisation
k = GPy.kern.rbf(d, lengthscale=2);      C = k.K(X,X); Z4 = np.random.multivariate_normal(mu,C,10)
k = GPy.kern.Matern32(d, lengthscale=2); C = k.K(X,X); Z5 = np.random.multivariate_normal(mu,C,10)
k = GPy.kern.Matern52(d, lengthscale=2); C = k.K(X,X); Z6 = np.random.multivariate_normal(mu,C,10)

# Export data to R for plotting
pandas.DataFrame(X).to_csv("../figures/ch4-example-gp-priorrealisations/ch4-example-gp-priorrealisations-x.csv");
pandas.DataFrame(Z1.transpose()).to_csv("../figures/ch4-example-gp-priorrealisations/ch4-example-gp-priorrealisations-z1.csv");
pandas.DataFrame(Z2.transpose()).to_csv("../figures/ch4-example-gp-priorrealisations/ch4-example-gp-priorrealisations-z2.csv");
pandas.DataFrame(Z3.transpose()).to_csv("../figures/ch4-example-gp-priorrealisations/ch4-example-gp-priorrealisations-z3.csv");
pandas.DataFrame(Z4.transpose()).to_csv("../figures/ch4-example-gp-priorrealisations/ch4-example-gp-priorrealisations-z4.csv");
pandas.DataFrame(Z5.transpose()).to_csv("../figures/ch4-example-gp-priorrealisations/ch4-example-gp-priorrealisations-z5.csv");
pandas.DataFrame(Z6.transpose()).to_csv("../figures/ch4-example-gp-priorrealisations/ch4-example-gp-priorrealisations-z6.csv");

########################################################################
# End of file
########################################################################