#!/usr/bin/env python
# Copyright (C) 2015 Andy Aschwanden

# nosetests orographic_precipy.py

import numpy as np
from scipy import integrate

# Set up constants in SI unuits
L = 30e3
nx =30e2
dx = 30e3

gamma = -6.5 * 1e-3

# Set up domain
x = np.linspace(0, L, nx+1)


def Ts(x):
    return TsL - gamma * z

def esat(T):
    '''
    Returns water vapor saturation pressure in Pa at temperature T [Celsius] according to
    August-Roche-Magnus formula
    '''
    
    e0 = 610.94
    a = 17.625
    b = 243.04
    
    return e0*np.exp((a*T)/(b+T))

def esat_test():
    
    T = np.array([-5, 0, 5], dtype=np.float64)

    assert np.max(np.fabs(esat(T) - np.array([ 421.9082468, 610.94, 871.55954949]))) < 1e-8


def gaussian_smoothing(t, x, dx):
    return exp(-((x-t)/dx)**2)

def minus_nabla_F():
    return alpha0 - alpha1*vbar*dzdx*esat

def integrand():
    return minus_nabla_F *gaussian_smoothing

def expint(x, dx):
    return quad(gaussian_smoothing, x, Inf, args=(x, dx))[0]

vec_expint = np.vectorize(expint)

