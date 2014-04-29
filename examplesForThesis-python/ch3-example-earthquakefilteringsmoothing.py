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
# State inference in the earthquake count model
# Example 3.2 in Section 3.3.2
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
sys              = stSystemEarthQuake()
sys.version      = "standard"
sys.par          = np.zeros((3,1))
sys.par[0]       = 0.88473937;
sys.par[1]       = 0.15055556;
sys.par[2]       = 17.64746228;
sys.T            = 114;

########################################################################
# Setup the parameters for the algorithm
########################################################################

par.fileprefix   = "earthquake"
par.nPars          = 1;
smc.nPart          = 1000;
smc.fixedLag       = 12;
smc.resampFactor   = 2;
smc.weightdist     = "poisson";

########################################################################
# Read the data
########################################################################

par.dataset = 0;
data.sample(sys,np.zeros(sys.T))
tmp = np.loadtxt("data/earthquake/earthquakes.txt"); data.y = tmp[:,1];

########################################################################
# Run the bPF to estimate the state and the predicted output
########################################################################

smc.bPF(data,sys,par)
plot(smc.xhatf)

out = hstack( ( smc.xhatf, sys.par[2] * np.exp(smc.xhatf) ) )
pandas.DataFrame(out).to_csv("ch3-example-earthquakefiltering.csv");

# Make some predictions on the future value of the data
xp = np.zeros((100,10));
xp[:,0] = smc.xhatf[sys.T-1];
for kk in range(0,100):
    for tt in range(0, 9):
        xp[kk,tt+1] = sys.f(xp[kk,tt], 0, 0, 0, tt) + sys.fn(0, 0, 0, tt) * np.random.randn(1);
    

# Export data to R for plotting
pandas.DataFrame(xp).to_csv("ch3-example-earthquakefiltering-pred.csv");

########################################################################
# End of file
########################################################################