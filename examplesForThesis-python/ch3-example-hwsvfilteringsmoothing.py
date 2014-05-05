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
# State inference in the Hull-White SV model
# Example 3.6 in Section 3.4.1
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
sys.par[0]       = 0.996;
sys.par[1]       = 0.129;
sys.par[2]       = 0.837;
sys.T            = 3567;

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
# Read the data
########################################################################

par.dataset = 0;
data.sample(sys,np.zeros(sys.T))
data.y = 100 * np.loadtxt("data/seOMXdata.csv",delimiter=",");

########################################################################
# Make grid and estimate ll, score and info
########################################################################

smc.flPS(data,sys,par)

# Plot the 95% CIs
plot( 1.96 * sys.par[2] * exp(0.5*smc.xhatf)) ; plot( 1.96 * sys.par[2] * exp(0.5*smc.xhats)); plot(data.y)

# Export data to R for plotting
out = hstack( (1.96 * sys.par[2] * exp(0.5*smc.xhatf), 1.96 * sys.par[2] * exp(0.5*smc.xhats), smc.xhatf, smc.xhats, np.matrix(data.y).transpose() ) )
pandas.DataFrame(out).to_csv("ch3-example-hwsvfilteringsmoothing.csv");

########################################################################
# End of file
########################################################################