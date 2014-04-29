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
# This file contains some helpers for the GPO algorithm
#
########################################################################

import numpy as np
import scipy as sp
import GPy

########################################################################
# EI rule
########################################################################

def EIeval(x,user_data):
    m = user_data[0];
    llmax = user_data[1];
    epsilon = user_data[2];
    
    # Find the predicted value in x
    Mup, Mus, up95, lo95 = m.predict( np.array(x) )
    
    # Calculate auxillary quantites
    s   = np.sqrt(Mus);
    yres  = Mup - llmax - epsilon;
    ynorm = ( yres / s) * ( s > 0 );
    
    # Compute the EI and negate it
    ei  = yres * sp.stats.norm.cdf(ynorm) + s * sp.stats.norm.pdf(ynorm);
    ei  = np.max((ei,0));
    return -ei, 0

########################################################################
# Evaluate the surrogate function to find its maxima
########################################################################

def MUeval(x,m):
    Mup, Mus, up95, lo95 = m.predict( np.array(x) )
    return -Mup, 0

########################################################################
# PI rule
########################################################################
    
def PIeval(x,user_data):
    m = user_data[0];
    llmax = user_data[1];
    epsilon = user_data[2];
    
    # Find the predicted value in x
    Mup, Mus, up95, lo95 = m.predict( np.array(x) )
    
    # Calculate auxillary quantites
    s   = np.sqrt(Mus);
    yres  = Mup - llmax - epsilon;
    ynorm = ( yres / s) * ( s > 0 );
    
    # Compute the EI and negate it
    pi  = sp.stats.norm.cdf(ynorm);
    return -pi, 0

########################################################################
# UCB rule
########################################################################

def UCBeval(x,user_data):
    m = user_data[0];
    llmax = user_data[1];
    epsilon = user_data[2];
    
    # Find the predicted value in x
    Mup, Mus, up95, lo95 = m.predict( np.array(x) )
    
    # Calculate auxillary quantites
    s   = np.sqrt(Mus);
    
    # Compute the EI and negate it
    ucb  = Mup + epsilon * s;
    return -ucb, 0

########################################################################
# End of file
########################################################################