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
# This file contains the Kalman filtering and smoothing algorithms
#
########################################################################

import numpy as np
import scipy as sp
from classes import *
from helpers import *

class kalmanFilter(object):
    ll      = [];
    xhatf   = [];
    xhatp   = [];    
    K       = [];
    Pf      = [];
    Pp      = [];
    filterType = "kalman";

    def filter(self,data,sys,par):
        A  = sys.par[0];
        #B  = sys.par[1];
        Q  = sys.par[1];
        C  = 1.0; #sys.par[3];
        #D  = sys.par[4];
        R  = sys.par[2];
        
        S       = np.zeros((sys.T,1));
        K       = np.zeros((sys.T,));
        xhatp   = np.zeros((sys.T+1,1));
        xhatf   = np.zeros((sys.T,1));
        yhatp   = np.zeros((sys.T,1));
        Pf      = np.zeros((sys.T,1));
        Pp      = np.zeros((sys.T+1,1));    
        ll = 0.0;
        
        Pp[0]     = par.Po;
        xhatp[0]  = par.xo;
        
        for tt in range(0, sys.T):
            # Calculate the Kalman Gain
            S[tt] = C * Pp[tt] * C + R;
            K[tt] = Pp[tt] * C / S[tt];
            
            # Compute the state estimate
            yhatp[tt]   = C * xhatp[tt];
            xhatf[tt]   = xhatp[tt] + K[tt] * ( data.y[tt] - yhatp[tt] );
            xhatp[tt+1] = A * xhatf[tt];
            
            # Update covariance
            Pf[tt]      = Pp[tt] - K[tt] * S[tt] * K[tt];
            Pp[tt+1]    = A * Pf[tt] * A + Q;
            
            # Estimate loglikelihood
            ll          -= 0.5 * np.log(2.0 * np.pi * S[tt]) + 0.5 * ( data.y[tt] - yhatp[tt] ) * ( data.y[tt] - yhatp[tt] ) / S[tt];
            
        self.ll    = ll;
        self.xhatf = xhatf;
        self.xhatp = xhatp;
        self.K     = K;
        self.Pp    = Pp;
        self.Pf    = Pf;

    def filtersd(self,data,sys,par):
        A  = sys.par[0];
        B  = 1.0;
        Q  = sys.par[1]**2;
        C  = 1.0; #sys.par[3];
        #D  = sys.par[4];
        R  = sys.par[2]**2;
        
        S       = np.zeros((sys.T,1));
        K       = np.zeros((sys.T,));
        xhatp   = np.zeros((sys.T+1,1));
        xhatf   = np.zeros((sys.T,1));
        yhatp   = np.zeros((sys.T,1));
        Pf      = np.zeros((sys.T,1));
        Pp      = np.zeros((sys.T+1,1));    
        ll = 0.0;
        
        Pp[0]     = par.Po;
        xhatp[0]  = par.xo;
        
        for tt in range(0, sys.T):
            # Calculate the Kalman Gain
            S[tt] = C * Pp[tt] * C + R;
            K[tt] = Pp[tt] * C / S[tt];
            
            # Compute the state estimate
            yhatp[tt]   = C * xhatp[tt];
            xhatf[tt]   = xhatp[tt] + K[tt] * ( data.y[tt] - yhatp[tt] );
            xhatp[tt+1] = A * xhatf[tt] + B * data.u[tt];
            
            # Update covariance
            Pf[tt]      = Pp[tt] - K[tt] * S[tt] * K[tt];
            Pp[tt+1]    = A * Pf[tt] * A + Q;
            
            # Estimate loglikelihood
            ll          -= 0.5 * np.log(2.0 * np.pi * S[tt]) + 0.5 * ( data.y[tt] - yhatp[tt] ) * ( data.y[tt] - yhatp[tt] ) / S[tt];
            
        self.ll    = ll;
        self.xhatf = xhatf;
        self.xhatp = xhatp;
        self.K     = K;
        self.Pp    = Pp;
        self.Pf    = Pf;
        
class kalmanSmoother(object):
    xhats   = [];
    score   = [];    
    info    = [];
    ll      = [];
    xhatf   = [];
    xhatp   = [];    
    K       = [];
    Pf      = [];
    Pp      = [];
        
    def RTS(self,data,sys,par):
        
        A  = sys.par[0];
        #B  = sys.par[1];
        Q  = sys.par[1];
        C  = 1.0; #sys.par[3];
        #D  = sys.par[4];
        R  = sys.par[2];
        
        # Run the Kalman filter
        S       = np.zeros((sys.T,1));
        K       = np.zeros((sys.T,1));
        xhatp   = np.zeros((sys.T+1,1));
        xhatf   = np.zeros((sys.T,1));
        yhatp   = np.zeros((sys.T,1));
        Pf      = np.zeros((sys.T,1));
        Pp      = np.zeros((sys.T+1,1));    
        ll = 0.0;
        
        Pp[0]     = par.Po;
        xhatp[0]  = par.xo;
        
        for tt in range(0, sys.T):
            # Calculate the Kalman Gain
            S[tt] = C * Pp[tt] * C + R;
            K[tt] = Pp[tt] * C / S[tt];
            
            # Compute the state estimate
            yhatp[tt]   = C * xhatp[tt];
            xhatf[tt]   = xhatp[tt] + K[tt] * ( data.y[tt] - yhatp[tt] );
            xhatp[tt+1] = A * xhatf[tt];
            
            # Update covariance
            Pf[tt]      = Pp[tt] - K[tt] * S[tt] * K[tt];
            Pp[tt+1]    = A * Pf[tt] * A + Q;
            
            # Estimate loglikelihood
            ll          -= 0.5 * np.log(2.0 * np.pi * S[tt]) + 0.5 * ( data.y[tt] - yhatp[tt] ) * ( data.y[tt] - yhatp[tt] ) / S[tt];
        
        # Run the Kalman RTS smoother
        J       = np.zeros((sys.T,1));
        M       = np.zeros((sys.T,1));
        xhats   = np.zeros((sys.T,1));
        Ps      = np.zeros((sys.T,1));
        
        Ps[sys.T-1]     = Pf[sys.T-1];
        xhats[sys.T-1]  = xhatf[sys.T-1];
        
        for tt in range((sys.T-2),0,-1):
            J[tt]       = Pf[tt] * A / Pp[tt+1]
            xhats[tt]   = xhatf[tt] + J[tt] * ( xhats[tt+1] - A * xhatf[tt] )
            Ps[tt]      = Pf[tt] + J[tt] * ( Ps[tt+1] - Pp[tt+1] ) * J[tt];
        
        # Calculate the M-matrix
        M[sys.T-1]  = ( 1 - K[sys.T-1] ) * A * Pf[sys.T-1];
        for tt in range((sys.T-2),0,-1):
            M[tt]   = Pf[tt] * J[tt-1] + J[tt-1] * ( M[tt+1] - A * Pf[tt] ) * J[tt-1];
        
        # Estimate the score
        s = np.zeros((3,1));
        
        for tt in range(1,sys.T):
            kappa = xhats[tt]   * data.y[tt];
            eta   = xhats[tt]   * xhats[tt]   + Ps[tt];
            eta1  = xhats[tt-1] * xhats[tt-1] + Ps[tt-1];
            psi   = xhats[tt-1] * xhats[tt]   + M[tt];
            
            s[0] += ( psi - A * eta ) / Q ;
            s[1] += 0.5 * ( A**2 * eta + eta1 - 2 * A * psi - Q ) / (Q**2);
            #s[3] += ( kappa - C * eta ) / R;
            #s[4] += 0.5 * ( C*C/(R*R) * eta + data.y[tt]**2/(R*R) - 1.0/R - 2.0/(R*R) * C * kappa);
        
        # Add the log-prior derivatives
        for nn in range(0,par.nPars):   
            s[nn]     = sys.dprior1(nn) + s[nn];
            
        self.ll    = ll;
        self.xhatf = xhatf;
        self.xhatp = xhatp;
        self.K     = K;
        self.Pp    = Pp;
        self.Pf    = Pf;
        self.xhats = xhats;
        self.score = s;
        
    def RTSsd(self,data,sys,par):
        
        A  = sys.par[0];
        #B  = sys.par[1];
        Q  = sys.par[1]**2;
        C  = 1.0; #sys.par[3];
        #D  = sys.par[4];
        R  = sys.par[2]**2;
        
        # Run the Kalman filter
        S       = np.zeros((sys.T,1));
        K       = np.zeros((sys.T,1));
        xhatp   = np.zeros((sys.T+1,1));
        xhatf   = np.zeros((sys.T,1));
        yhatp   = np.zeros((sys.T,1));
        Pf      = np.zeros((sys.T,1));
        Pp      = np.zeros((sys.T+1,1));    
        ll = 0.0;
        
        Pp[0]     = par.Po;
        xhatp[0]  = par.xo;
        
        for tt in range(0, sys.T):
            # Calculate the Kalman Gain
            S[tt] = C * Pp[tt] * C + R;
            K[tt] = Pp[tt] * C / S[tt];
            
            # Compute the state estimate
            yhatp[tt]   = C * xhatp[tt];
            xhatf[tt]   = xhatp[tt] + K[tt] * ( data.y[tt] - yhatp[tt] );
            xhatp[tt+1] = A * xhatf[tt];
            
            # Update covariance
            Pf[tt]      = Pp[tt] - K[tt] * S[tt] * K[tt];
            Pp[tt+1]    = A * Pf[tt] * A + Q;
            
            # Estimate loglikelihood
            ll          -= 0.5 * np.log(2.0 * np.pi * S[tt]) + 0.5 * ( data.y[tt] - yhatp[tt] ) * ( data.y[tt] - yhatp[tt] ) / S[tt];
        
        # Run the Kalman RTS smoother
        J       = np.zeros((sys.T,1));
        M       = np.zeros((sys.T,1));
        xhats   = np.zeros((sys.T,1));
        Ps      = np.zeros((sys.T,1));
        
        Ps[sys.T-1]     = Pf[sys.T-1];
        xhats[sys.T-1]  = xhatf[sys.T-1];
        
        for tt in range((sys.T-2),0,-1):
            J[tt]       = Pf[tt] * A / Pp[tt+1]
            xhats[tt]   = xhatf[tt] + J[tt] * ( xhats[tt+1] - A * xhatf[tt] )
            Ps[tt]      = Pf[tt] + J[tt] * ( Ps[tt+1] - Pp[tt+1] ) * J[tt];
        
        # Calculate the M-matrix
        M[sys.T-1]  = ( 1 - K[sys.T-1] ) * A * Pf[sys.T-1];
        for tt in range((sys.T-2),0,-1):
            M[tt]   = Pf[tt] * J[tt-1] + J[tt-1] * ( M[tt+1] - A * Pf[tt] ) * J[tt-1];
        
        # Estimate the score
        s = np.zeros((3,1));
        
        for tt in range(1,sys.T):
            kappa = xhats[tt]   * data.y[tt];
            eta   = xhats[tt]   * xhats[tt]   + Ps[tt];
            eta1  = xhats[tt-1] * xhats[tt-1] + Ps[tt-1];
            psi   = xhats[tt-1] * xhats[tt]   + M[tt];
            
            s[0] += ( psi - A * eta ) / Q ;
            s[1] += 0.5 * ( A**2 * eta + eta1 - 2 * A * psi - Q ) / (Q**2);
            #s[3] += ( kappa - C * eta ) / R;
            #s[4] += 0.5 * ( C*C/(R*R) * eta + data.y[tt]**2/(R*R) - 1.0/R - 2.0/(R*R) * C * kappa);
        
        # Add the log-prior derivatives
        for nn in range(0,par.nPars):   
            s[nn]     = sys.dprior1(nn) + s[nn];
            
        self.ll    = ll;
        self.xhatf = xhatf;
        self.xhatp = xhatp;
        self.K     = K;
        self.Pp    = Pp;
        self.Pf    = Pf;
        self.xhats = xhats;
        self.score = s;    

########################################################################
# End of file
########################################################################