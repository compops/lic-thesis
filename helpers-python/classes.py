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
# This file contains the different models and data generation
#
########################################################################

import numpy as np

########################################################################
# Class definitions
########################################################################

class stParameters(object):
    nPars                 = 1;
    nMCMC                 = 1000;
    nBurnIn               = 500;
    stepSize              = 0.1;
    
    normLimit             = 0.01;
    dataset               = 0.0;
    fileprefix            = [];
    nProgressReport       = 100;
    verboseSampler        = 0;
    zeroscore             = 0;
    writeOutPriorWarnings = 0;

class stSystemLGSS(object):
    ########################################################################
    # LGSS model
    ########################################################################
    #
    # x(t+1) = par[0] * x(t) + par[1] * v(t),    v(t) ~ N(0,1)
    # y(t)   =          x(t) + par[2] * e(t),    e(t) ~ N(0,1)
    #
    ########################################################################

    nPar       = 3;
    par        = np.zeros(3);
    T          = 0
    xo         = 0
    so         = 0;
    model      = "Linear Gaussian system"
    supportsFA = 1.0;
    scale      = 1.0;
    version    = "standard"
    dh         = 0.0;

    def storeParameters(self,newParm,sys,par):
        self.par = np.zeros(self.nPar);

        for kk in range(0,par.nPars):
            self.par[kk] = np.array(newParm[kk], copy=True)

        for kk in range(par.nPars,self.nPar):
            self.par[kk] = sys.par[kk];

    def returnParameters(self,par):
        out = np.zeros(par.nPars);

        for kk in range(0,par.nPars):
            out[kk]  = self.par[kk];

        return(out);

    def transform(self):
        if (self.version == "tanhexp"):
            self.par[0] = np.tanh( self.par[0] );
            self.par[1] = np.exp ( self.par[1] ) * self.scale;
        else:
            self.par[1] = self.par[1] * self.scale;

    def invTransform(self):
        if (self.version == "tanhexp"):
            self.par[0] = np.arctanh( self.par[0] );
            self.par[1] = np.log    ( self.par[1] / self.scale );
        else:
            self.par[1] = self.par[1] / self.scale;

    def h(self, xt, st, ut, tt):
        return 0.0;
    
    def f(self, xt, ut, st, yt, tt):
        return self.par[0] * xt + ut;
    
    def fn(self, xt, st, yt, tt):
        return self.par[1];        
    
    def g(self, xt, ut, st, tt):
        return xt;

    def gn(self, xt, st, tt):
        return self.par[2];

    def fa(self, xt, ytt, ut, st, tt):
        delta = self.par[1]**(-2) + self.par[2]**(-2); delta = 1.0 / delta;
        return delta * ( ytt * self.par[2]**(-2) + self.par[1]**(-2) * self.par[0] * xt );

    def fna(self, xt, ytt, st, tt):
        delta = self.par[1]**(-2) + self.par[2]**(-2); delta = 1.0 / delta;
        return np.sqrt(delta);
    
    def ga(self, xt, ut, st, tt):
        return self.par[0] * xt;

    def gna(self, xt, st, tt):
        return np.sqrt( self.par[1]**2 + self.par[2]**2 );   
    
    def Dparm(self, xtt, xt, yt, st, at, ut, par):
        
        nOut = len(xtt);
        gradient = np.zeros(( nOut, par.nPars ));
        Q1 = self.par[1]**(-1);        
        Q2 = self.par[1]**(-2);
        Q3 = self.par[1]**(-3);        
        R1 = self.par[2]**(-1);
        R3 = self.par[2]**(-3);
        px = xtt - self.par[0] * xt - ut;
        py = yt - xt;

        if (self.version == "tanhexp"):
            for v1 in range(0,par.nPars):
                if v1 == 0:
                    gradient[:,v1] = xt * Q2 * px * ( 1.0 - self.par[0]**2 );
                elif v1 == 1:
                    gradient[:,v1] = ( Q2 * px**2 - 1.0 ) * self.scale;
                elif v1 == 2:
                    gradient[:,v1] = R3 * py**2 - R1;
                else:
                    gradient[:,v1] = 0.0;
        else:
            for v1 in range(0,par.nPars):
                if v1 == 0:
                    gradient[:,v1] = xt * Q2 * px;
                elif v1 == 1:
                    gradient[:,v1] = ( Q3 * px**2 - Q1 ) * self.scale;
                elif v1 == 2:
                    gradient[:,v1] = R3 * py**2 - R1;
                else:
                    gradient[:,v1] = 0.0;            
        return(gradient);


    def DDparm(self, xtt, xt, yt, st, at, ut, par):
        
        nOut = len(xtt);
        hessian = np.zeros( (nOut, par.nPars,par.nPars) );
        Q1 = self.par[1]**(-1);
        Q2 = self.par[1]**(-2);
        Q3 = self.par[1]**(-3);
        Q4 = self.par[1]**(-4);
        R2 = self.par[2]**(-2);
        R4 = self.par[2]**(-4);
        px = xtt - self.par[0] * xt;
        py = yt - xt;

        if (self.version == "tanhexp"):
            for v1 in range(0,par.nPars):
                for v2 in range(0,par.nPars):
                    if ( (v1 == 0) & (v2 == 0) ):
                        hessian[:,v1,v2] = - xt**2 * Q2 * ( 1.0 - self.par[0]**2 )**2 - 2.0 * self.par[0] * Q2 * xt * px * ( 1.0 - self.par[0]**2 )
    
                    elif ( (v1 == 1) & (v2 == 1) ):
                        hessian[:,v1,v2] = - 2.0 * Q2 * px**2 * self.scale**2;
    
                    elif ( ( (v1 == 1) & (v2 == 0) ) | ( (v1 == 0) & (v2 == 1) ) ):
                        hessian[:,v1,v2] = - 2.0 * xt * Q2 * px * ( 1.0 - self.par[0] ) * self.scale;
    
                    elif ( (v1 == 2) & (v2 == 2) ):
                        hessian[:,v1,v2] = R2 - 3.0 * R4 * py**2
    
                    else:
                        hessian[:,v1,v2] = 0.0;
                    
        else:
            for v1 in range(0,par.nPars):
                for v2 in range(0,par.nPars):
                    if ( (v1 == 0) & (v2 == 0) ):
                        hessian[:,v1,v2] = - xt**2 * Q2;
    
                    elif ( (v1 == 1) & (v2 == 1) ):
                        hessian[:,v1,v2] = ( Q2 - 3.0 * Q4 * px**2 - Q1 ) * self.scale**2;
    
                    elif ( ( (v1 == 1) & (v2 == 0) ) | ( (v1 == 0) & (v2 == 1) ) ):
                        hessian[:,v1,v2] = - 2.0 * xt * Q3 * px * self.scale;
    
                    elif ( (v1 == 2) & (v2 == 2) ):
                        hessian[:,v1,v2] = R2 - 3.0 * R4 * py**2
    
                    else:
                        hessian[:,v1,v2] = 0.0;            

        return(hessian);

    def priorUniform(self):
        return 1.0;

    def prior(self):
        return(0.0);

    def dprior1(self, v1):
        return(0.0);


    def ddprior1(self, v1, v2):
        return(0.0);

    def Jacobian( self ):
        if (self.version == "tanhexp"):
            return np.log( 1.0 - self.par[0]**2 ) + np.log( self.par[1] );
        else:
            return 0.0;

class stSystemGARCH(object):
    ########################################################################
    # GARCH(1,1) in noise model
    ########################################################################
    #
    # s(t+1) = par[0] + par[1] * x(t)**2 + par[2] * s(t)
    # x(t+1) = s(t+1) * v(t),                               v(t) ~ N(0,1)
    # y(t)   = x(t) + par[4] * e(t),                        e(t) ~ N(0,1)
    #
    ########################################################################

    nPar       = 4;
    par        = np.zeros(4);
    T          = 0
    xo         = 0
    so         = 0
    model      = "GARCH(1,1) in noise"
    supportsFA = 1.0;

    def storeParameters(self,newParm,sys,par):
        self.par = np.zeros(self.nPar);

        for kk in range(0,par.nPars):
            self.par[kk] = np.array(newParm[kk], copy=True)

        for kk in range(par.nPars,self.nPar):
            self.par[kk] = sys.par[kk];

    def returnParameters(self,par):
        out = np.zeros(par.nPars);

        for kk in range(0,par.nPars):
            out[kk]  = self.par[kk];

        return(out);

    def transform(self):
        self.par[0] = self.par[0]
        self.par[1] = self.par[1]

    def invTransform(self):
        self.par[0] = self.par[0]
        self.par[1] = self.par[1]

    def h(self, xt, ut, st, tt):
        return self.par[0] + self.par[1] * xt**2 + self.par[2] * st;
    
    def f(self, xt, ut, st, yt, tt):
        return 0;
    
    def fn(self, xt, st, yt, tt):
        return np.sqrt( st );   
    
    def g(self, xt, ut, st, tt):
        return xt;

    def gn(self, xt, st, tt):
        return self.par[3];

    def fa(self, xt, ytt, ut, st, tt):
        delta = st**(-1) + self.par[3]**(-2);
        return delta**(-1) * ytt * self.par[3]**(-2);

    def fna(self, xt, ytt, st, tt):
        delta = st**(-1) + self.par[3]**(-2);
        return np.sqrt( delta**(-1) );
    
    def ga(self, xt, ut, st, tt):
        return 0;

    def gna(self, xt, st, tt):
        return np.sqrt( self.par[3]**2 + st );   
    
    def Dparm(self, xtt, xt, yt, st, at, ut, par):
        
        nOut = len(xtt);
        gradient = np.zeros(( nOut, par.nPars ));      
        h1 = ( self.par[0] + self.par[1] * xt**2 + self.par[2] * st )**(-1);
        h2 = ( self.par[0] + self.par[1] * xt**2 + self.par[2] * st )**(-2);
        R1 = self.par[3]**(-1);
        R3 = self.par[3]**(-3);
        py = yt - xt;
        
        self.rD1[:,0] = 1     +  self.par[2] * self.rD1[at,0];
        self.rD1[:,1] = xt**2 +  self.par[2] * self.rD1[at,1];
        self.rD1[:,2] = st    +  self.par[2] * self.rD1[at,2];

        for v1 in range(0,par.nPars):
            if v1 == 0:
                gradient[:,v1] = 0.5 * ( h2 * xtt**2 - h1 ) * self.rD1[:,0];
            elif v1 == 1:
                gradient[:,v1] = 0.5 * ( h2 * xtt**2 - h1 ) * self.rD1[:,1];
            elif v1 == 2:
                gradient[:,v1] = 0.5 * ( h2 * xtt**2 - h1 ) * self.rD1[:,2];
            elif v1 == 3:
                gradient[:,v1] = R3 * py**2 - R1;
            else:
                gradient[:,v1] = 0.0;  
        
        return(gradient);


    def DDparm(self, xtt, xt, yt, st, at, ut, par):
        
        nOut = len(xtt);
        hessian = np.zeros( (nOut, par.nPars,par.nPars) );
        h1 = ( self.par[0] + self.par[1] * xt**2 + self.par[2] * st )**(-1);
        h2 = ( self.par[0] + self.par[1] * xt**2 + self.par[2] * st )**(-2);
        h3 = ( self.par[0] + self.par[1] * xt**2 + self.par[2] * st )**(-3);
        R2 = self.par[3]**(-2);
        R4 = self.par[3]**(-4);
        py = yt - xt;

        # Recursions for the second derivatives of log f wrt. to (alpha, beta, gamma) and gamma (the others are zero)
        self.rDD2[:,0] =       self.rDD1[at,0] +               self.rDD2[at,0];
        self.rDD2[:,1] =       self.rDD1[at,1] +               self.rDD2[at,1];
        self.rDD2[:,2] = 2.0 * self.rDD1[at,2] + self.par[2] * self.rDD2[at,2]
        
        # Recursions for the first derivatives of log f wrt. to (alpha, beta, gamma)
        self.rDD1[:,0] = 1     +  self.par[2] * self.rDD1[at,0];
        self.rDD1[:,1] = xt**2 +  self.par[2] * self.rDD1[at,1];
        self.rDD1[:,2] = st    +  self.par[2] * self.rDD1[at,2];
        

        for v1 in range(0,par.nPars):
            for v2 in range(0,par.nPars):
                if ( (v1 == 0) & (v2 == 0) ):
                    hessian[:,v1,v2] = 0.5 * ( h2 - 2.0 * h3 * xtt**2 ) * self.rDD1[:,0]**2 ;

                elif ( (v1 == 1) & (v2 == 1) ):
                    hessian[:,v1,v2] = 0.5 * ( h2 - 2.0 * h3 * xtt**2 ) * self.rDD1[:,1]**2 ;
                
                elif ( (v1 == 2) & (v2 == 2) ):
                    hessian[:,v1,v2] = 0.5 * ( h2 - 2.0 * h3 * xtt**2 ) * self.rDD1[:,2]**2 + 0.5 * ( h2 * xtt**2 - h1 ) * self.rDD2[:,2];
                    
                elif ( (v1 == 3) & (v2 == 3) ):
                    hessian[:,v1,v2] = R2 - 3.0 * R4 * py**2;
                    
                elif ( ( (v1 == 1) & (v2 == 0) ) | ( (v1 == 0) & (v2 == 1) ) ):
                    hessian[:,v1,v2] = 0.5 * ( h2 - 2.0 * h3 * xtt**2 ) * self.rDD1[:,0]*self.rDD1[:,1];

                elif ( ( (v1 == 2) & (v2 == 0) ) | ( (v1 == 0) & (v2 == 2) ) ):
                    hessian[:,v1,v2] = 0.5 * ( h2 - 2.0 * h3 * xtt**2 ) * self.rDD1[:,0]*self.rDD1[:,2] + 0.5 * ( h2 * xtt**2 - h1 ) * self.rDD2[:,0];

                elif ( ( (v1 == 1) & (v2 == 2) ) | ( (v1 == 2) & (v2 == 1) ) ):
                    hessian[:,v1,v2] = 0.5 * ( h2 - 2.0 * h3 * xtt**2 ) * self.rDD1[:,1]*self.rDD1[:,2] + 0.5 * ( h2 * xtt**2 - h1 ) * self.rDD2[:,1];
                
                else:
                    hessian[:,v1,v2] = 0.0;            
        
        return(hessian);

    def priorUniform(self):
        out = 1.0;

        if( self.par[0] < 0.0 ):
            out = 0.0;
                
        if( self.par[1] < 0.0 ):
            out = 0.0;

        if( self.par[2] < 0.0 ):
            out = 0.0;

        if( ( self.par[1] + self.par[2] ) > 1.0 ):
            out = 0.0;
            
        if( self.par[3] < 0.0 ):
            out = 0.0;
        
        return( out );
    
    def prior(self):
        return(0.0);

    def dprior1(self, v1):
        return(0.0);

    def ddprior1(self, v1, v2):
        return(0.0);

    def Jacobian( self ):
        return 0.0;

class stSystemHW(object):
    ########################################################################
    # HW model
    ########################################################################
    #
    # x(t+1) = par[0] * x(t) + par[1] * v(t),                v(t) ~ N(0,1)
    # y(t)   =                 par[2] exp(-x(t)/2)* e(t),    e(t) ~ N(0,1)
    #
    ########################################################################

    nPar       = 3;
    par        = np.zeros(3);
    T          = 0
    xo         = 0
    so         = 0;
    model      = "Hull-White Stochastic Volatility model"
    supportsFA = 0.0;
    scale      = 1.0;
    version    = "standard"
    dh         = 0.0;

    def storeParameters(self,newParm,sys,par):
        self.par = np.zeros(self.nPar);

        for kk in range(0,par.nPars):
            self.par[kk] = np.array(newParm[kk], copy=True)

        for kk in range(par.nPars,self.nPar):
            self.par[kk] = sys.par[kk];

    def returnParameters(self,par):
        out = np.zeros(par.nPars);

        for kk in range(0,par.nPars):
            out[kk]  = self.par[kk];

        return(out);

    def transform(self):
        if (self.version == "tanhexp"):
            self.par[0] = np.tanh( self.par[0] );
            self.par[1] = np.exp ( self.par[1] ) * self.scale;
        else:
            self.par[1] = self.par[1] * self.scale;

    def invTransform(self):
        if (self.version == "tanhexp"):
            self.par[0] = np.arctanh( self.par[0] );
            self.par[1] = np.log    ( self.par[1] / self.scale );
        else:
            self.par[1] = self.par[1] / self.scale;

    def h(self, xt, st, ut, tt):
        return 0.0;
    
    def f(self, xt, ut, st, yt, tt):
        return self.par[0] * xt;
    
    def fn(self, xt, st, yt, tt):
        return self.par[1];        
    
    def g(self, xt, ut, st, tt):
        return 0;

    def gn(self, xt, st, tt):
        return self.par[2] * np.exp( 0.5 * xt);

    def fa(self, xt, ytt, ut, st, tt):
        return 0;

    def fna(self, xt, ytt, st, tt):
        return 0;
    
    def ga(self, xt, ut, st, tt):
        return 0;

    def gna(self, xt, st, tt):
        return 0;
    
    def Dparm(self, xtt, xt, yt, st, at, ut, par):
        
        nOut = len(xtt);
        gradient = np.zeros(( nOut, par.nPars ));
        Q1 = self.par[1]**(-1);        
        Q2 = self.par[1]**(-2);
        Q3 = self.par[1]**(-3);        
        R1 = self.par[2]**(-1);
        R3 = self.par[2]**(-3);
        px = xtt - self.par[0] * xt;
        py = yt;

        if (self.version == "tanhexp"):
            for v1 in range(0,par.nPars):
                if v1 == 0:
                    gradient[:,v1] = xt * Q2 * px * ( 1.0 - self.par[0]**2 );
                elif v1 == 1:
                    gradient[:,v1] = ( Q2 * px**2 - 1.0 ) * self.scale;
                elif v1 == 2:
                    gradient[:,v1] = ( R3 * py**2 - R1 ) * np.exp( 0.5 * xt );
                else:
                    gradient[:,v1] = 0.0;
        else:
            for v1 in range(0,par.nPars):
                if v1 == 0:
                    gradient[:,v1] = xt * Q2 * px;
                elif v1 == 1:
                    gradient[:,v1] = ( Q3 * px**2 - Q1 ) * self.scale;
                elif v1 == 2:
                    gradient[:,v1] = ( R3 * py**2 - R1 ) * np.exp( 0.5 * xt );
                else:
                    gradient[:,v1] = 0.0;            
        return(gradient);


    def DDparm(self, xtt, xt, yt, st, at, ut, par):
        
        nOut = len(xtt);
        hessian = np.zeros( (nOut, par.nPars,par.nPars) );
        Q1 = self.par[1]**(-1);
        Q2 = self.par[1]**(-2);
        Q3 = self.par[1]**(-3);
        Q4 = self.par[1]**(-4);
        R2 = self.par[2]**(-2);
        R4 = self.par[2]**(-4);
        px = xtt - self.par[0] * xt;
        py = yt;

        if (self.version == "tanhexp"):
            for v1 in range(0,par.nPars):
                for v2 in range(0,par.nPars):
                    if ( (v1 == 0) & (v2 == 0) ):
                        hessian[:,v1,v2] = - xt**2 * Q2 * ( 1.0 - self.par[0]**2 )**2 - 2.0 * self.par[0] * Q2 * xt * px * ( 1.0 - self.par[0]**2 )
    
                    elif ( (v1 == 1) & (v2 == 1) ):
                        hessian[:,v1,v2] = - 2.0 * Q2 * px**2 * self.scale**2;
    
                    elif ( ( (v1 == 1) & (v2 == 0) ) | ( (v1 == 0) & (v2 == 1) ) ):
                        hessian[:,v1,v2] = - 2.0 * xt * Q2 * px * ( 1.0 - self.par[0] ) * self.scale;
    
                    elif ( (v1 == 2) & (v2 == 2) ):
                        hessian[:,v1,v2] = ( R2 - 3.0 * R4 * py**2 ) * np.exp( xt )
    
                    else:
                        hessian[:,v1,v2] = 0.0;
                    
        else:
            for v1 in range(0,par.nPars):
                for v2 in range(0,par.nPars):
                    if ( (v1 == 0) & (v2 == 0) ):
                        hessian[:,v1,v2] = - xt**2 * Q2;
    
                    elif ( (v1 == 1) & (v2 == 1) ):
                        hessian[:,v1,v2] = ( Q2 - 3.0 * Q4 * px**2 - Q1 ) * self.scale**2;
    
                    elif ( ( (v1 == 1) & (v2 == 0) ) | ( (v1 == 0) & (v2 == 1) ) ):
                        hessian[:,v1,v2] = - 2.0 * xt * Q3 * px * self.scale;
    
                    elif ( (v1 == 2) & (v2 == 2) ):
                        hessian[:,v1,v2] = ( R2 - 3.0 * R4 * py**2 ) * np.exp( xt )
    
                    else:
                        hessian[:,v1,v2] = 0.0;            

        return(hessian);

    def priorUniform(self):
        return 1.0;

    def prior(self):
        return(0.0);

    def dprior1(self, v1):
        return(0.0);


    def ddprior1(self, v1, v2):
        return(0.0);

    def Jacobian( self ):
        if (self.version == "tanhexp"):
            return np.log( 1.0 - self.par[0]**2 ) + np.log( self.par[1] );
        else:
            return 0.0;

class stSystemEarthQuake(object):
    ########################################################################
    # Earthquake model
    ########################################################################

    nPar       = 3;
    par        = np.zeros(3);
    T          = 0
    xo         = 0
    so         = 0;
    model      = "Earthquake model"
    supportsFA = 0.0;
    scale      = 1.0;
    version    = "standard"
    dh         = 0.0;

    def storeParameters(self,newParm,sys,par):
        self.par = np.zeros(self.nPar);

        for kk in range(0,par.nPars):
            self.par[kk] = np.array(newParm[kk], copy=True)

        for kk in range(par.nPars,self.nPar):
            self.par[kk] = sys.par[kk];

    def returnParameters(self,par):
        out = np.zeros(par.nPars);

        for kk in range(0,par.nPars):
            out[kk]  = self.par[kk];

        return(out);

    def transform(self):
        return 0.0;

    def invTransform(self):
        return 0.0;

    def h(self, xt, st, ut, tt):
        return 0.0;
    
    def f(self, xt, ut, st, yt, tt):
        return self.par[0] * xt;
    
    def fn(self, xt, st, yt, tt):
        return self.par[1];        
    
    def g(self, xt, ut, st, tt):
        return self.par[2] * np.exp( xt );

    def gn(self, xt, st, tt):
        return 0;

    def fa(self, xt, ytt, ut, st, tt):
        return 0;
    
    def fna(self, xt, ytt, st, tt):
        return 0;
    
    def ga(self, xt, ut, st, tt):
        return 0;

    def gna(self, xt, st, tt):
        return 0;
    
    def Dparm(self, xtt, xt, yt, st, at, ut, par):
        return 0;

    def DDparm(self, xtt, xt, yt, st, at, ut, par):
        return 0;

    def priorUniform(self):
        return 1.0;

    def prior(self):
        return(0.0);

    def dprior1(self, v1):
        return(0.0);


    def ddprior1(self, v1, v2):
        return(0.0);

    def Jacobian( self ):
        return 0.0;

class stSystemWeather(object):
    ########################################################################
    # Earthquake model
    ########################################################################

    nPar       = 7;
    par        = np.zeros(7);
    T          = 0
    xo         = 0
    so         = 0;
    model      = "weather model"
    supportsFA = 0.0;
    scale      = 1.0;
    version    = "standard"
    dh         = 0.0;

    def storeParameters(self,newParm,sys,par):
        self.par = np.zeros(self.nPar);

        for kk in range(0,par.nPars):
            self.par[kk] = np.array(newParm[kk], copy=True)

        for kk in range(par.nPars,self.nPar):
            self.par[kk] = sys.par[kk];

    def returnParameters(self,par):
        out = np.zeros(par.nPars);

        for kk in range(0,par.nPars):
            out[kk]  = self.par[kk];

        return(out);

    def transform(self):
        return 0.0;

    def invTransform(self):
        return 0.0;

    def h(self, xt, st, ut, tt):
        return self.par[3] * np.cos( 2 * tt * np.pi / 365) + self.par[4] * np.cos( 4 * tt * np.pi / 365) + self.par[5] * np.sin( 2 * tt * np.pi / 365) + self.par[6] * np.sin( 4 * tt * np.pi / 365);
    
    def f(self, xt, ut, st, yt, tt):
        return self.par[0] * xt;
    
    def fn(self, xt, st, yt, tt):
        return self.par[1];        
    
    def g(self, xt, ut, st, tt):
        ll = np.exp( self.par[2] + xt + st);
        return ll / ( 1 + ll );

    def gn(self, xt, st, tt):
        return 0;

    def fa(self, xt, ytt, ut, st, tt):
        return 0;
    
    def fna(self, xt, ytt, st, tt):
        return 0;
    
    def ga(self, xt, ut, st, tt):
        return 0;

    def gna(self, xt, st, tt):
        return 0;
    
    def Dparm(self, xtt, xt, yt, st, at, ut, par):
        return 0;

    def DDparm(self, xtt, xt, yt, st, at, ut, par):
        return 0;

    def priorUniform(self):
        return 1.0;

    def prior(self):
        return(0.0);

    def dprior1(self, v1):
        return(0.0);


    def ddprior1(self, v1, v2):
        return(0.0);

    def Jacobian( self ):
        return 0.0;

########################################################################
# Data class for generating samples
########################################################################

class stData(object):
    x = []
    u = []
    y = []
    e = []
    v = []
    T = []
    model = []

    # the class "constructor"  - It's actually an initializer
    def __init__(self):
        model = "empty"

    def sample(self, sys, u):
        s = np.zeros((sys.T+1,1));
        x = np.zeros((sys.T+1,1));
        y = np.zeros((sys.T,1));
        v = np.random.randn(sys.T);
        e = np.random.randn(sys.T);

        x[0] = sys.xo;
        s[0] = sys.so;       

        for tt in range(0, sys.T):
            y[tt]   = sys.g(x[tt], u[tt], s[tt],tt) + sys.gn(x[tt], s[tt], tt) * e[tt];
            s[tt+1] = sys.h(x[tt], u[tt], s[tt], tt);
            x[tt+1] = sys.f(x[tt], u[tt], s[tt], y[tt], tt) + sys.fn(x[tt], s[tt], y[tt], tt) * v[tt];

        self.x = x[0:sys.T]
        self.s = s[0:sys.T]
        self.y = y;
        self.u = u;
        self.v = v;
        self.e = e;
        self.T = sys.T;
        self.model = sys.model;

########################################################################
# End of file
########################################################################