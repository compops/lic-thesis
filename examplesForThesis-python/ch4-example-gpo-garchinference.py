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
# GPO for ML inference in the GARCH(1,1) model
# Example 4.8 in Section 4.4.3 
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
sys              = stSystemGARCH()
sys.version      = "standard"
sys.par          = np.zeros((4,1))
sys.par[0]       = 0.10;
sys.par[1]       = 0.80;
sys.par[2]       = 0.05;
sys.par[3]       = 0.30;
sys.T            = 250;

thSys            = stSystemGARCH()
thSys.T          = sys.T;
thSys.version    = sys.version;

# Generate some data
data.sample(sys, np.zeros(sys.T));

########################################################################
# Setup the parameters for the algorithm
########################################################################

par.fileprefix     = "garch"
par.nPars          = 1;
par.epsilonEI      = 0.01;
par.epsilonUCB     = 0.01;
par.preIter        = 2;
par.maxIter        = 10;

smc.nPart          = 100;
smc.resamplingType = "systematic";     # multinomial or systematic
smc.filterType     = "fullyadapted";   # kalman or bootstrap or fullyadapted

########################################################################
# Preliminaries
########################################################################

# Make a grid to evaluate the EI on
l      = np.array([0.01], dtype=np.float64)
u      = np.array([0.25], dtype=np.float64)

# Preallocate the information needed
th     = np.zeros((par.maxIter+1,par.nPars))
EI     = np.zeros((par.maxIter+1,1))
mumax  = np.zeros((par.maxIter,1))
pSamp  = np.zeros((par.maxIter,1))
thhat  = np.zeros((par.maxIter,par.nPars))
llmax  = np.zeros((par.maxIter,1))
grid  = np.arange(0.01,0.26,0.01);
gridM = np.zeros((len(grid),par.maxIter));
gridL = np.zeros((len(grid),par.maxIter));
gridU = np.zeros((len(grid),par.maxIter));
gridA = np.zeros((len(grid),par.maxIter));

#-----------------------------------------------------------------------
# Evaluate the log-likelihood in some random points (many) to get the hyperparameters
#-----------------------------------------------------------------------
tmpth = np.zeros((30,1));
tmpll = np.zeros((30,1));
for kk in range(0,30):
    tmpth[kk]     =  l + (u-l) * np.random.random(par.nPars)
    
    thSys.storeParameters((tmpth[kk], sys.par[1], sys.par[2], sys.par[3]),sys,par);
    smc.faPF(data, thSys, par)
    tmpll[kk] = smc.ll;
    print kk;

# Fit the hyperparameters of the kernel
kernel = GPy.kern.Matern52(input_dim=par.nPars) + GPy.kern.bias(input_dim=par.nPars)
m = GPy.models.GPRegression(tmpth,tmpll,kernel)
m.optimize()
m.optimize_restarts(num_restarts = 10)  

########################################################################
# GPO using the EI rule
########################################################################

# Make a preliminary run to estimate the log-likelihood
for kk in range(0,par.preIter):
    th[kk,:]     =  l + (u-l) * np.random.random(par.nPars)
    
    thSys.storeParameters((th[kk,:], sys.par[1], sys.par[2], sys.par[3]),sys,par);
    smc.faPF(data, thSys, par)
    pSamp[kk] = smc.ll;
    print kk;

#-----------------------------------------------------------------------
# Main GPO-loop
#-----------------------------------------------------------------------
for kk in range(par.preIter-1,par.maxIter):

    # Sample the log-likelihood
    thP = (th[kk,:], sys.par[1], sys.par[2], sys.par[3]);
    thSys.storeParameters(thP,sys,par);
    smc.faPF(data, thSys, par)
    pSamp[kk] = smc.ll;
    
    # Fit the GP
    thPart      = th[0:kk+1,:];
    pSampPart   = pSamp[0:kk+1];
    m = GPy.models.GPRegression(thPart,pSampPart,kernel)
    
    # Find the maximum expected value
    Mup, ys2, up95, lo95 = m.predict( np.array(thPart) )
    mumax[kk] = np.max(Mup);
    
    # Compute the next point in which to sample the posterior
    x, fmin, ierror = solve(EIeval,l,u,user_data=(m,mumax[kk],par.epsilonEI),maxf=1000,maxT=1000);
    
    # Set the new point and save the estimate of the EI
    th[kk+1,:] = x;
    EI[kk+1] = fmin;
    
    tmp2, tmp3, ierror = solve(MUeval,l,u,user_data=m,algmethod=1);
    thhat[kk,:] = tmp2;
    llmax[kk,:] = tmp3;
    
    user_data = ( m,mumax[kk],par.epsilonEI)
    Mup, ys2, up95, lo95 = m.predict( np.array(grid.reshape((len(grid),1))) )
    
    for ii in range(0,len(grid)):    
        gridA[ii,kk], tmp = EIeval( np.array(grid[ii]), user_data )
    
    gridM[:,kk] = Mup[:,0];
    gridL[:,kk] = lo95[:,0];
    gridU[:,kk] = up95[:,0];
    
    print kk;

# Export data to R for plotting
out = hstack( ( gridM, gridL, gridU, gridA  ) )
pandas.DataFrame(out).to_csv("ch4-example-gpo-garchinference-ei.csv");

# Export data to R for plotting
out = hstack( ( th.transpose(), pSamp.transpose()  ) )
pandas.DataFrame(out).to_csv("ch4-example-gpo-garchinference-ei2.csv");

########################################################################
# GPO using the PI rule
########################################################################

# Preallocate the information needed
th     = np.zeros((par.maxIter+1,par.nPars))
EI     = np.zeros((par.maxIter+1,1))
mumax  = np.zeros((par.maxIter,1))
pSamp  = np.zeros((par.maxIter,1))
thhat  = np.zeros((par.maxIter,par.nPars))
llmax  = np.zeros((par.maxIter,1))
grid  = np.arange(0.01,0.26,0.01);
gridM = np.zeros((len(grid),par.maxIter));
gridL = np.zeros((len(grid),par.maxIter));
gridU = np.zeros((len(grid),par.maxIter));
gridA = np.zeros((len(grid),par.maxIter));

# Main GPO-loop
for kk in range(par.preIter-1,par.maxIter):

    # Sample the log-likelihood
    thP = (th[kk,:], sys.par[1], sys.par[2], sys.par[3]);
    thSys.storeParameters(thP,sys,par);
    smc.faPF(data, thSys, par)
    pSamp[kk] = smc.ll;
    
    # Fit the GP
    thPart      = th[0:kk+1,:];
    pSampPart   = pSamp[0:kk+1];
    m = GPy.models.GPRegression(thPart,pSampPart,kernel)
    
    # Find the maximum expected value
    Mup, ys2, up95, lo95 = m.predict( np.array(thPart) )
    mumax[kk] = np.max(Mup);
    
    # Compute the next point in which to sample the posterior
    x, fmin, ierror = solve(EIeval,l,u,user_data=(m,mumax[kk],par.epsilonEI),maxf=1000,maxT=1000);
    
    # Set the new point and save the estimate of the EI
    th[kk+1,:] = x;
    EI[kk+1] = fmin;
    
    tmp2, tmp3, ierror = solve(MUeval,l,u,user_data=m,algmethod=1);
    thhat[kk,:] = tmp2;
    llmax[kk,:] = tmp3;
    
    user_data = ( m,mumax[kk],par.epsilonEI)
    Mup, ys2, up95, lo95 = m.predict( np.array(grid.reshape((len(grid),1))) )
    
    for ii in range(0,len(grid)):    
        gridA[ii,kk], tmp = PIeval( np.array(grid[ii]), user_data )
    
    gridM[:,kk] = Mup[:,0];
    gridL[:,kk] = lo95[:,0];
    gridU[:,kk] = up95[:,0];
    
    print kk;

# Export data to R for plotting
out = hstack( ( gridM, gridL, gridU, gridA  ) )
pandas.DataFrame(out).to_csv("ch4-example-gpo-garchinference-pi.csv");

# Export data to R for plotting
out = hstack( ( th.transpose(), pSamp.transpose()  ) )
pandas.DataFrame(out).to_csv("ch4-example-gpo-garchinference-pi2.csv");

########################################################################
# GPO using the UCB rule
########################################################################

# Preallocate the information needed
th     = np.zeros((par.maxIter+1,par.nPars))
EI     = np.zeros((par.maxIter+1,1))
mumax  = np.zeros((par.maxIter,1))
pSamp  = np.zeros((par.maxIter,1))
thhat  = np.zeros((par.maxIter,par.nPars))
llmax  = np.zeros((par.maxIter,1))
grid  = np.arange(0.01,0.26,0.01);
gridM = np.zeros((len(grid),par.maxIter));
gridL = np.zeros((len(grid),par.maxIter));
gridU = np.zeros((len(grid),par.maxIter));
gridA = np.zeros((len(grid),par.maxIter));

# Main GPO-loop
for kk in range(par.preIter-1,par.maxIter):

    # Sample the log-likelihood
    thP = (th[kk,:], sys.par[1], sys.par[2], sys.par[3]);
    thSys.storeParameters(thP,sys,par);
    smc.faPF(data, thSys, par)
    pSamp[kk] = smc.ll;
    
    # Fit the GP
    thPart      = th[0:kk+1,:];
    pSampPart   = pSamp[0:kk+1];
    m = GPy.models.GPRegression(thPart,pSampPart,kernel)
    
    # Find the maximum expected value
    Mup, ys2, up95, lo95 = m.predict( np.array(thPart) )
    mumax[kk] = np.max(Mup);
    
    # Compute the next point in which to sample the posterior
    x, fmin, ierror = solve(EIeval,l,u,user_data=(m,mumax[kk],par.epsilonEI),maxf=1000,maxT=1000);
    
    # Set the new point and save the estimate of the EI
    th[kk+1,:] = x;
    EI[kk+1] = fmin;
    
    tmp2, tmp3, ierror = solve(MUeval,l,u,user_data=m,algmethod=1);
    thhat[kk,:] = tmp2;
    llmax[kk,:] = tmp3;
    
    user_data = ( m,mumax[kk],1.96)
    Mup, ys2, up95, lo95 = m.predict( np.array(grid.reshape((len(grid),1))) )
    
    for ii in range(0,len(grid)):    
        gridA[ii,kk], tmp = UCBeval( np.array(grid[ii]), user_data )
    
    gridM[:,kk] = Mup[:,0];
    gridL[:,kk] = lo95[:,0];
    gridU[:,kk] = up95[:,0];
    
    print kk;

# Export data to R for plotting
out = hstack( ( gridM, gridL, gridU, gridA  ) )
pandas.DataFrame(out).to_csv("ch4-example-gpo-garchinference-ucb.csv");

# Export data to R for plotting
out = hstack( ( th.transpose(), pSamp.transpose()  ) )
pandas.DataFrame(out).to_csv("ch4-example-gpo-garchinference-ucb2.csv");

########################################################################
# End of file
########################################################################