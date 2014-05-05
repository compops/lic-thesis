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
# The IS algorithm for parameter inference
# Example 3.1 in Section 3.2 
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
sys.par[0]       = 0.50;
sys.par[1]       = 1.00;
sys.par[2]       = 0.10;
sys.T            = 3567;

thSys            = stSystemHW()
thSys.T          = sys.T;
thSys.par        = np.zeros((3,1));

########################################################################
# Setup the parameters for the algorithm
########################################################################

par.fileprefix   = "hwsv"
par.nPars        = 3;
par.xo           = 0;
par.Po           = 0.001;

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
# Read the data
########################################################################

par.dataset = 0;
data.sample(sys,np.zeros(sys.T))
data.y = 100 * np.loadtxt("data/seOMXdata.csv",delimiter=",");

########################################################################
# Importance sampling
########################################################################

# No. samples from the target
N = 2000;

# Allocate vectors
th     = np.zeros((N,3))
ll     = np.zeros(N)

# For each sample
for kk in range(0,N):
    # propose
    th[kk,0] = np.random.uniform(0.50,1.00,1);
    th[kk,1] = np.random.gamma(2,0.1,1);
    th[kk,2] = np.random.gamma(7,0.1,1);
        
    # compute target
    thSys.par[0] = th[kk,0];
    thSys.par[1] = th[kk,1];
    thSys.par[2] = th[kk,2];
    
    smc.flPS(data,thSys,par);
    ll[kk]    = smc.ll;

# Export data to R for plotting
out = vstack((th[:,0],th[:,1],th[:,2],ll))
pandas.DataFrame(out).to_csv("ch3-example-importancesampling-hwsv.csv");

########################################################################
# End of file
########################################################################