#!/usr/bin/env python
# Copyright (C) 2015 Andy Aschwanden

# nosetests orographic_precipy.py

import numpy as np
from scipy import integrate

# Set up constants in SI units
L = 30e3  # m
nx =30e2  #
dx = 30e3  # m

gamma = -6.5 * 1e-3  # degC / 1000 m

vbar = 0.5  # m/s

TsL = 5  # degC

# Set up domain
x = np.linspace(-L/2, L/2, nx+1)

def gaussian_bump(x, a=1/(2*np.pi), b=0, c=1):
    
    return a*np.exp((-(x-b)**2)/(2*c**2))

def gaussian_bump_test():
    
    sigma = np.sqrt(0.5)
    mu = -2
    a = 1/(sigma*np.sqrt(2*np.pi))
    b = mu
    c = sigma
    x = np.linspace(-5, 2, 8)
    assert  np.max(np.fabs(gaussian_bump(x, a, b, c) - np.array([  6.96265260e-05, 1.03334927e-02,
                                                                2.07553749e-01, 5.64189584e-01,
                                                                2.07553749e-01,   1.03334927e-02,
                                                                6.96265260e-05,   6.34911734e-08]))) < 1e-9

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

    assert  np.max(np.fabs(esat(T) - np.array([ 421.9082468, 610.94, 871.55954949]))) < 1e-8


def gaussian_smoothing(t, x, dx):
    return exp(-((x-t)/dx)**2)

def alphas():

    return 1. / esat(T), 110  # m/yr/

def minus_nabla_F():
    
    return alpha0 - alpha1*vbar*dzdx*esat

def integrand():
    return minus_nabla_F *gaussian_smoothing

def expint(x, dx):
    return quad(gaussian_smoothing, x, Inf, args=(x, dx))[0]

vec_expint = np.vectorize(expint)

def hillslope_erosion(x, t):
    from scipy.special import erfc
    D = 1
    return 0.5 * z0 * erfc(x / (2 * np.sqrt(D * t)))

import pylab as plt

z = gaussian_bump(x, 2e3, 0, L/2)
plt.plot(x, z)

