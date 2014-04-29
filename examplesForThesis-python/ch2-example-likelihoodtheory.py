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
# Log-likelihood, score and information matrix for an LGSS model
# Example 2.6 in Section 2.3 
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
ks               = kalmanSmoother();
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
# Read the data from file
########################################################################

par.dataset = 0;
file = 'data/' + str(par.fileprefix) + 'T' + str(sys.T) +'/' + str(par.fileprefix) + 'DataT' + str(sys.T) + str(par.dataset) + '.csv'
tmp = np.loadtxt(file,delimiter=",")
data.x = tmp[:,0]; data.u = tmp[:,1]; data.y = tmp[:,2];

########################################################################
# Make grid and estimate ll, score and info
########################################################################

# Make grid
xx = arange(0.00,1.00,0.01);

# Allocate arrays
ll     = np.zeros(len(xx))
score  = np.zeros(len(xx))
info   = np.zeros(len(xx))
ddatay = np.zeros((100,sys.T));
sscore = np.zeros((100,1));

# Loop over grid to estimate the log-likelihood and score
for ii in range(0,len(xx)):
    
    # Set the parameters and run the RTS smoother
    sys.par[0] = xx[ii];
    ks.RTSsd(data,sys,par)
    
    ll[ii]    = ks.ll;
    score[ii] = ks.score[0];
    print((ii,len(xx)))

# Generate 100 data sets from the model
sys.par[0] = 0.5;
for jj in range(0,100):
    data.sample(sys,np.zeros(sys.T))
    ddatay[jj,:] = data.y[:,0];

# Compute the score function for each data set and grid point
for ii in range(0,len(xx)):
        sys.par[0] = xx[ii];
        for jj in range(0,100):
            data.y[:,0] = ddatay[jj,:];
            ks.RTSsd(data,sys,par);
            sscore[jj] = smc.score[0];
        
        # Compute the expected information as the sample covariance
        info[ii] = np.var(sscore);
        print((ii,len(xx)))
        
# Plot some graphs
subplot(3,1,1); plot(xx,ll)
subplot(3,1,2); plot(xx,score)
subplot(3,1,3); plot(xx,info)

# Export data to R for plotting
out = vstack((xx,ll,score,info)).transpose()
pandas.DataFrame(out).to_csv("ch2-example-likelihoodtheory.csv");

########################################################################
# End of file
########################################################################