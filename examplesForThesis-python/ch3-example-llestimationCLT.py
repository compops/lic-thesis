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
# Bias and varaiance of the log-likelihood estimate
# Example 3.4 in Section 3.3.4
#
########################################################################

from smc import *
from kf import *
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
kf               = kalmanFilter();

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

smc.nPart          = 10;
smc.resamplingType = "systematic";           # multinomial or systematic
smc.filterType     = "fullyadapted";         # kalman or bootstrap or fullyadapted
smc.smootherType   = "fixedlag";             # kalman or filtersmoother or fixedlag or ffbsm (not implemented)
smc.flVersion      = "full";                 # filtersmoother or neglectcross or full
smc.fixedLag       = 12;
smc.onlydiagInfo   = 0;
smc.makeInfoPSD    = 1;
smc.resampFactor   = 2;

########################################################################
# Generate some data
########################################################################

par.dataset = 0;
data.sample(sys,np.zeros(sys.T))

########################################################################
# Make grid and estimate ll, score and info
########################################################################

llPF     = np.zeros(1000)

# Compute 1000 estimates of the log-likelihood using the faPF on the same data
for ii in range(0,1000):
    smc.faPF(data,sys,par)     
    llPF[ii]    = smc.ll;    

# Compute the true log-likelihood with the Kalman filter
kf.filtersd(data,sys,par)
llTrue = kf.ll;

# Export data to R for plotting
out = vstack((llPF-llTrue))
pandas.DataFrame(out).to_csv("ch3-example-llestimationCLT.csv");

########################################################################
# End of file
########################################################################
