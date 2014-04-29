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
# Score and information matrix in the Hull-White SV model
# Example 3.7 in Section 3.4.2
#
########################################################################

from smc import *
from classes import *
from helpers import *
import pandas
import numpy as np

########################################################################
# Arrange the data structures
########################################################################
data             = stData();
smc              = smcSampler();
par              = stParameters();

########################################################################
# Setup the system
########################################################################
sys              = stSystemHW()
sys.version      = "standard"
sys.par          = np.zeros((3,1))
sys.par[0]       = 0.95;
sys.par[1]       = 0.16;
sys.par[2]       = 0.70;
sys.T            = 250;

########################################################################
# Setup the parameters for the algorithm
########################################################################

par.fileprefix   = "hw"
par.nPars          = 1;

smc.nPart          = 5000;
smc.resamplingType = "systematic";           # multinomial or systematic
smc.filterType     = "bootstrap";            # kalman or bootstrap or fullyadapted
smc.smootherType   = "fixedlag";             # kalman or filtersmoother or fixedlag or ffbsm (not implemented)
smc.flVersion      = "full";                 # filtersmoother or neglectcross or full
smc.fixedLag       = 12;
smc.onlydiagInfo   = 0;
smc.makeInfoPSD    = 1;
smc.resampFactor   = 2;

########################################################################
# Generate data from the model
########################################################################

par.dataset = 0;
data.sample(sys,np.zeros(sys.T))

########################################################################
# Make grid and estimate ll, score and info
########################################################################

# Make the grid
xx = arange(0.70,1.00,0.005);

# Allocate vectors
ll     = np.zeros(len(xx))
score  = np.zeros(len(xx))
info   = np.zeros(len(xx))
ng     = np.zeros(len(xx))

# Run the FL smoother for each grid point
for ii in range(0,len(xx)):
    sys.par[0] = xx[ii];
    smc.flPS(data,sys,par)
     
    ll[ii]    = smc.ll;
    score[ii] = smc.score;
    info[ii]  = smc.infom;
    
    print((ii,len(xx)))

# Plot the score and observed information matrix
subplot(2,1,1); plot(xx,score)
subplot(2,1,2); plot(xx,info)

# Export data to R for plotting
out = vstack((xx,ll,score,info)).transpose()
pandas.DataFrame(out).to_csv("ch3-example-hwsvllscoreinfo.csv");

########################################################################
# End of file
########################################################################
