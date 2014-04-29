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
# Parameter inference in the LGSS model
# Example 4.1 in Section 4.2
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
sys.T            = 250;

thSys            = stSystemLGSS()
thSys.T          = sys.T;
thSys.version    = sys.version

########################################################################
# Setup the parameters for the algorithm
########################################################################

par.fileprefix   = "lgss"
par.nPars          = 2;
par.nMCMC          = 10000;
par.nBurnIn        = 500;
par.normLimit      = 0.1;
par.verboseSampler = 0;

########################################################################
# Select the initial point for the sampler
########################################################################

par.initPar      = np.zeros((3,1))
par.initPar[0]   = 0.1;
par.initPar[1]   = 2.0;
par.initPar[2]   = 0.1;
par.Po           = 0.01;
par.xo           = 0;

########################################################################
# Read the data
########################################################################

par.dataset = 0;
file = 'data/' + str(par.fileprefix) + 'T' + str(sys.T) +'/' + str(par.fileprefix) + 'DataT' + str(sys.T) + str(par.dataset) + '.csv'
tmp = np.loadtxt(file,delimiter=",")
data.x = tmp[:,0]; data.u = tmp[:,1]; data.y = tmp[:,2];

########################################################################
# Run the sampler
########################################################################
par.stepSize = (0.10, 0.10, 0.10);
mh.runSampler(kf, data, sys, thSys, par, "PMH0");

subplot(2,1,1); plot(mh.th[:,0])
subplot(2,1,2); plot(mh.th[:,1])

np.mean(mh.accept)

# Export data to R for plotting
out = hstack((mh.th,mh.accept))
pandas.DataFrame(out).to_csv("ch4-example-lgssinference.csv");

########################################################################
# End of file
########################################################################