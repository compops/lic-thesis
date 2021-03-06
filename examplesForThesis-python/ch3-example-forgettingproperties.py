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
# Mixing property in the LGSS model
# Example 3.5 in Section 3.4.1
#
########################################################################

from kf import *
from pmh import *
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
mh               = stPMH();

########################################################################
# Setup the system
########################################################################
sys              = stSystemLGSS()
sys.version      = "standard"
sys.par          = np.zeros((3,1))
sys.par[0]       = 0.50;
sys.par[1]       = 1.00;
sys.par[2]       = 0.10;
sys.T            = 15;

thSys            = stSystemLGSS()
thSys.T          = sys.T;
thSys.version    = sys.version

########################################################################
# Setup the parameters for the algorithm
########################################################################

par.fileprefix   = "lgss"
par.nPars          = 1;
par.normLimit      = 0.1;
par.verboseSampler = 0;

par.Po           = 0.01;
par.xo           = 0;

########################################################################
# Run the sampler
########################################################################

nRuns = 20;
xx = np.zeros((nRuns,sys.T));

# Mixing at phi=0.20
sys.par[0]       = 0.20;
for ii in range(0,nRuns):
    
    # Generate initial state
    sys.xo = np.random.normal(0,20);
    
    # Generate data and save the state
    data.sample(sys,np.zeros(sys.T))
    xx[ii,:] = (data.x).transpose();

# Export data to R for plotting
plot(xx.transpose())
pandas.DataFrame(xx.transpose()).to_csv("ch3-example-forgettingproperties1.csv");

sys.par[0]       = 0.50;
for ii in range(0,nRuns):
    
    # Generate initial state
    sys.xo = np.random.normal(0,20);
    
    # Generate data and save the state
    data.sample(sys,np.zeros(sys.T))
    xx[ii,:] = (data.x).transpose();

# Export data to R for plotting
plot(xx.transpose())
pandas.DataFrame(xx.transpose()).to_csv("ch3-example-forgettingproperties2.csv");

sys.par[0]       = 0.80;
for ii in range(0,nRuns):
    
    # Generate initial state
    sys.xo = np.random.normal(0,20);
    
    # Generate data and save the state
    data.sample(sys,np.zeros(sys.T))
    xx[ii,:] = (data.x).transpose();

# Export data to R for plotting
plot(xx.transpose())
pandas.DataFrame(xx.transpose()).to_csv("ch3-example-forgettingproperties3.csv");

########################################################################
# End of file
########################################################################
