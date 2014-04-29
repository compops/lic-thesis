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
# Path degeneracy in the GARCH(1,1) model
# Example 3.3 in Section 3.3.2
#
########################################################################

from smc import *
from pmh import *
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
sys              = stSystemGARCH()
sys.par          = np.zeros((4,1))
sys.par[0]       = 0.10;
sys.par[1]       = 0.80;
sys.par[2]       = 0.05;
sys.par[3]       = 0.30;
sys.T            = 21;

########################################################################
# Setup the parameters for the algorithm
########################################################################

par.fileprefix      = "garch"
par.nPars           = 4;

smc.nPart          = 10;
smc.resamplingType = "systematic";     # multinomial or systematic
smc.filterType     = "bootstrap";   # kalman or bootstrap or fullyadapted
smc.smootherType   = "fixedlag";       # kalman or filtersmoother or fixedlag or ffbsm (not implemented)
smc.flVersion      = "full";           # multinomial or systematic# filtersmoother or neglectcross or full
smc.fixedLag       = 16;
smc.onlydiagInfo   = 0;
smc.makeInfoPSD    = 1;
par.nProgressReport= 5;

########################################################################
# Read the data from file
########################################################################

par.dataset = 0;
file = 'data/garch/simulateddata.csv'
tmp = np.loadtxt(file,delimiter=",")
data.y = tmp[0:sys.T]; data.u = np.zeros(sys.T);

########################################################################
# Run the filters and export the particles and ancestors to R for plotting
########################################################################

smc.nPart          = 10;
smc.filterType     = "bootstrap";
smc.flPS(data,sys,par)
pandas.DataFrame(smc.p).to_csv("ch3-example-garchpathdegenercy-bpf10-p.csv");
pandas.DataFrame(smc.a).to_csv("ch3-example-garchpathdegenercy-bpf10-a.csv");

smc.nPart          = 20;
smc.filterType     = "bootstrap";
smc.flPS(data,sys,par)
pandas.DataFrame(smc.p).to_csv("ch3-example-garchpathdegenercy-bpf20-p.csv");
pandas.DataFrame(smc.a).to_csv("ch3-example-garchpathdegenercy-bpf20-a.csv");

smc.nPart          = 10;
smc.filterType     = "fullyadapted";
smc.flPS(data,sys,par)
pandas.DataFrame(smc.p).to_csv("ch3-example-garchpathdegenercy-fapf-p.csv");
pandas.DataFrame(smc.a).to_csv("ch3-example-garchpathdegenercy-fapf-a.csv");

########################################################################
# End of file
########################################################################