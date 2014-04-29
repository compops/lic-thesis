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
# PMH0 for parameter inference in the GARCH(1,1) model
# Example 4.4 in Section 4.3
#
########################################################################

from smc import *
from pmh import *
from classes import *
from helpers import *
import pandas
import numpy as np
import os

########################################################################
# Arrange the data structures
########################################################################
data             = stData();
smc              = smcSampler();
par              = stParameters();
pmh0             = stPMH();
pmh1             = stPMH();
pmh2             = stPMH();
pmh0r            = stPMH();
pmh1r            = stPMH();
pmh2r            = stPMH();

########################################################################
# Setup the system
########################################################################
sys              = stSystemGARCH()
sys.par          = np.zeros((4,1))
sys.par[0]       = 0.1;
sys.par[1]       = 0.2;
sys.par[2]       = 0.2;
sys.par[3]       = 1.0;
sys.T            = 3567;

thSys            = stSystemGARCH()
thSys.T          = sys.T;

########################################################################
# Setup the parameters for the algorithm
########################################################################

par.fileprefix      = "garch"
par.nPars           = 4;

par.nMCMC          = 20000;
par.nBurnIn        = 0;
par.normLimit      = 0.1;
par.verboseSampler = 0;

smc.nPart          = 3000;
smc.resamplingType = "systematic";     # multinomial or systematic
smc.filterType     = "fullyadapted";   # kalman or bootstrap or fullyadapted
smc.smootherType   = "fixedlag";       # kalman or filtersmoother or fixedlag or ffbsm (not implemented)
smc.flVersion      = "full";           # multinomial or systematic# filtersmoother or neglectcross or full
smc.fixedLag       = 16;
smc.onlydiagInfo   = 0;
smc.makeInfoPSD    = 1;
par.nProgressReport= 50;

########################################################################
# Select the initial point for the sampler
########################################################################

par.initPar      = np.zeros((4,1))
par.initPar[0]   = sys.par[0];
par.initPar[1]   = sys.par[1];
par.initPar[2]   = sys.par[2];
par.initPar[3]   = sys.par[3];

########################################################################
# Read the data
########################################################################

par.dataset = 0;
file = 'data/seOMXdata.csv'; tmp = np.loadtxt(file,delimiter=",")
data.y = 100 * tmp; data.u = np.zeros(sys.T);

########################################################################
# Run the samplers
########################################################################
par.stepSize = (0.005, 0.005, 0.005);
pmh0.runSampler(smc, data, sys, thSys, par,"PMH0");

# Export data to R for plotting
pmh0.writeToFile("results/garch/omx/pmh0-run.csv",par);

########################################################################
# End of file
########################################################################
        
