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
# Input design in the LGSS model using GPO
# Example 4.10 in Section 4.4.3
#
########################################################################

from smc import *
from classes import *
from helpers import *
import numpy as np
from DIRECT import solve
import GPy
from gpohelpers import *

########################################################################
# Arrange the data structures
########################################################################
data             = stData();
smc              = smcSampler();
par              = stParameters();

########################################################################
# Setup the system
########################################################################
sys              = stSystemLGSS()
sys.version      = "standard"
sys.par          = np.zeros((3,1))
sys.par[0]       = 0.50;
sys.par[1]       = 0.10;
sys.par[2]       = 0.10;
sys.T            = 250;

########################################################################
# Setup the parameters for the algorithm
########################################################################

par.fileprefix     = "lgss"
par.nPars          = 2;
par.epsilonEI      = 0.01;
par.preIter        = 10;
par.maxIter        = 50;

smc.nPart          = 100;
smc.resamplingType = "systematic";     # multinomial or systematic
smc.filterType     = "fullyadapted";   # kalman or bootstrap or fullyadapted
smc.smootherType   = "fixedlag";       # kalman or filtersmoother or fixedlag or ffbsm (not implemented)
smc.flVersion      = "full";           # filtersmoother or neglectcross or full
smc.fixedLag       = 12;
smc.onlydiagInfo   = 0;
smc.makeInfoPSD    = 0;

########################################################################
# Preliminaries
########################################################################

# No. samples to estimate the information matrix from
nRep        = 100;

# Make a grid to evaluate the EI on
l      = np.array([0.01, 0.01], dtype=np.float64)
u      = np.array([0.99, 0.99], dtype=np.float64)

# Preallocate the information needed
th     = np.zeros((par.maxIter+1,par.nPars))
EI     = np.zeros((par.maxIter+1,1))
mumax  = np.zeros((par.maxIter,1))
pSamp  = np.zeros((par.maxIter,1))
score  = np.zeros((nRep,par.nPars))

# Compute the expected information matrix in some random points
for kk in range(0,par.preIter):
    th[kk,:]     = np.random.random(2)
    
    # Compute an input using the current parameters
    inputdata = 1 * ( th[kk,0] < np.random.uniform(size=sys.T) ) * np.sign( np.random.uniform(size=sys.T) - th[kk,1] )
    
    # For each data realisation
    for ii in range(0,nRep):
            
        # Create a data realisation
        data.sample(sys, inputdata )
        
        # Compute the score vector
        smc.flPS(data, sys, par)
        score[ii,] = smc.score;
        
    # Estimate the expected information matrix using the sample covariance
    pSamp[kk] = np.log( np.linalg.det( np.cov(score.transpose()) ) )
    print kk;

# Fit the hyperparameters of the GP kernel
kernel = GPy.kern.Matern52(input_dim=par.nPars, ARD=True) + GPy.kern.bias(input_dim=par.nPars)
m = GPy.models.GPRegression(th[0:par.preIter,:],pSamp[0:par.preIter],kernel)
m.optimize()
m.optimize_restarts(num_restarts = 10)

########################################################################
# Main GPO-loop
########################################################################

for kk in range(par.preIter-1,par.maxIter):

    # Sample the expected information matrix
    thP = th[kk,0];
    
    # Compute an input using the current parameters
    inputdata = 1 * ( th[kk,0] < np.random.uniform(size=sys.T) ) * np.sign( np.random.uniform(size=sys.T) - th[kk,1] )
    
    # For each data realisation
    for ii in range(0,nRep):
        
        # Create a data realisation
        data.sample(sys, inputdata )
        
        # Compute the score vector
        smc.flPS(data, sys, par)
        score[ii,] = smc.score;
    
    # Estimate the expected information matrix using the sample covariance
    pSamp[kk] = np.log( np.linalg.det( np.cov(score.transpose()) ) )
    
    # Fit the GP
    thPart      = th[0:kk+1,:];
    pSampPart   = pSamp[0:kk+1];
    m = GPy.models.GPRegression(thPart,pSampPart,kernel)
        
    if (( np.remainder(kk,10) == 0) & ( kk > par.preIter) ) :
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
    
    print kk;

# Estimate the argument that maximising the expected information matrix
muhat, llmax, ierror = solve(MUeval,l,u,user_data=m);

# Evaluate the surrogate function over a grid
xx = arange(0.00,1.01,0.01)
yy = arange(0.00,1.01,0.01)
mout =  np.zeros((len(xx)*len(yy)))
kk=0;
for ii in range(0,len(xx)):
    for jj in range(0,len(yy)):
        tmp, ys2, up95, lo95 = m.predict( np.array( (xx[ii],yy[jj]) ) )
        mout[kk] = tmp;
        kk=kk+1;
    
# Export data to R for plotting
import pandas
pandas.DataFrame(mout).to_csv("ch4-example-lgss-inputdesign.csv");
pandas.DataFrame(th).to_csv("ch4-example-lgss-inputdesign-th.csv");
pandas.DataFrame(pSamp).to_csv("ch4-example-lgss-inputdesign-pSamp.csv");

########################################################################
# End of file
########################################################################