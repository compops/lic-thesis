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
# GPO for ML inference in the earthquake count model
# Example 4.9 in Section 4.4.3
#
########################################################################

from smc import *
from classes import *
from helpers import *
import numpy as np
from DIRECT import solve
import GPy
from gpohelpers import *
import pandas

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
sys.par[0]       = 0.89;
sys.par[1]       = 0.14;
sys.par[2]       = 17.81;
sys.T            = 114;

thSys            = stSystemEarthQuake()
thSys.T          = sys.T;
thSys.version    = sys.version;

# Read the data
data.sample(sys, np.zeros(sys.T));
tmp = np.loadtxt("data/earthquake/earthquakes.txt"); data.y = tmp[:,1];

########################################################################
# Setup the parameters for the algorithm
########################################################################

par.fileprefix     = "earth"
par.nPars          = 3;
par.epsilonEI      = 0.01;          # Epsilon in the EI rule
par.preIter        = 50;            # No. initial samples from the GPO alg.
par.maxIter        = 200;           # No. iterations of the GPO alg.

smc.nPart          = 1000;
smc.resamplingType = "systematic";  # multinomial or systematic
smc.filterType     = "bootstrap";   # kalman or bootstrap or fullyadapted
smc.weightdist     = "poisson";

########################################################################
# Main GPO-loop
########################################################################

# Make a grid to evaluate the EI on
l      = np.array([0.01,0.01,10.0], dtype=np.float64)
u      = np.array([0.90,1.00,20.0], dtype=np.float64)

# Preallocate the information needed
th     = np.zeros((par.maxIter+1,par.nPars))
EI     = np.zeros((par.maxIter+1,1))
mumax  = np.zeros((par.maxIter,1))
pSamp  = np.zeros((par.maxIter,1))
thhat  = np.zeros((par.maxIter,3))
llmax  = np.zeros((par.maxIter,1))

# Make a preliminary run to estimate expected information is some random points
for kk in range(0,par.preIter):
    th[kk,:]     = l + (u-l) * np.random.random(par.nPars)
    
    thSys.storeParameters(th[kk,:],sys,par);
    smc.bPF(data, thSys, par)
    pSamp[kk] = smc.ll;
    
    print kk;

# Compute the kernel hyperparameters using the random samples
kernel = GPy.kern.Matern52(input_dim=par.nPars, ARD=True) + GPy.kern.bias(input_dim=par.nPars)
m = GPy.models.GPRegression(th[0:par.preIter,:],pSamp[0:par.preIter],kernel)
m.optimize()
m.optimize_restarts(num_restarts = 10)
        
# Main GPO-loop
for kk in range(par.preIter-1,par.maxIter):

    # Sample the log-likelihood
    thP = th[kk,:];
    thSys.storeParameters(th[kk,:],sys,par);
    smc.bPF(data, thSys, par)
    pSamp[kk] = smc.ll;
    
    # Fit the GP
    thPart      = th[0:kk+1,:];
    pSampPart   = pSamp[0:kk+1];
    m = GPy.models.GPRegression(thPart,pSampPart,kernel)

    if ( np.remainder(kk,10) == 0):
        m.optimize()
        m.optimize_restarts(num_restarts = 10)
    
    # Find the maximum expected value
    Mup, ys2, up95, lo95 = m.predict( np.array(thPart) )
    mumax[kk] = np.max(Mup);
    
    # Compute the next point in which to sample the posterior
    x, fmin, ierror = solve(EIeval,l,u,user_data=(m,mumax[kk],par.epsilonEI),maxf=1000,maxT=1000);
    
    # Set the new point and save the estimate of the EI
    th[kk+1,:] = x;
    EI[kk+1] = fmin;
    
    # Compute the current ML parameter estimate
    tmp2, tmp3, ierror = solve(MUeval,l,u,user_data=m,algmethod=1);
    thhat[kk,:] = tmp2;
    llmax[kk,:] = tmp3;
    
    print kk;

# Export data to R for plotting
out = hstack( ( thhat, llmax ) )
pandas.DataFrame(out).to_csv("ch4-example-earthinference.csv");

########################################################################
# End of file
########################################################################