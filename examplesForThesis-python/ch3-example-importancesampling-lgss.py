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
# IS for Bayesian parameter inference in the LGSS model
# Example 3.1 in Section 3.2 
#
########################################################################

from kf import *
from classes import *
from helpers import *
import pandas
import numpy as np

########################################################################
# Arrange the data structures
########################################################################
data             = stData();
kf               = kalmanFilter();
par              = stParameters();

########################################################################
# Setup the system
########################################################################
sys              = stSystemLGSS()
sys.version      = "standard"
sys.par          = np.zeros((3,1))
sys.par[0]       = 0.50;
sys.par[1]       = 1.00;
sys.par[2]       = 0.10;
sys.T            = 250;

########################################################################
# Setup the parameters for the algorithm
########################################################################

par.fileprefix   = "lgss"
par.nPars        = 1;
par.xo           = 0;
par.Po           = 0.001;

########################################################################
# Read the data
########################################################################

par.dataset = 0;
file = 'data/' + str(par.fileprefix) + 'T' + str(sys.T) +'/' + str(par.fileprefix) + 'DataT' + str(sys.T) + str(par.dataset) + '.csv'
tmp = np.loadtxt(file,delimiter=",")
data.x = tmp[:,0]; data.u = tmp[:,1]; data.y = tmp[:,2];

########################################################################
# Importance sampling
########################################################################

# No. samples from the target
N = 500;

# Allocate vectors
th     = np.zeros((N,2))
ll     = np.zeros(N)

# For each sample
for kk in range(0,N):
    # propose
    th[kk,0] = np.random.uniform(-0.99,0.99,1);
    th[kk,1] = np.random.gamma(1,1,1);
        
    # compute target
    sys.par[0] = th[kk,0];
    sys.par[1] = th[kk,1];
    kf.filtersd(data,sys,par);
    ll[kk] = kf.ll;

# Compute the normalised weights
normfactor = sum(exp(ll));
w = exp(ll) / normfactor;

# Compute the parameter estimates
sum( th[:,0] * w )
sum( th[:,1] * w )

# Export data to R for plotting
out = vstack((th[:,0],th[:,1],w))
pandas.DataFrame(out).to_csv("ch3-example-importancesampling-lgss.csv");

########################################################################
# End of file
########################################################################
