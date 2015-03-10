#!/usr/bin/env python
# Copyright (C) 2015 Andy Aschwanden

import numpy as np
from scipy import integrate

# Set up constants in SI unuits
L = 30e3
nx =30e2
dx = 30e3

gamma = -6.5 * 1e-3
e0 =

# Set up domain
x = np.linspace(0, L, nx+1)


def Ts(x):
    return TsL - gamma * z

def esat(Ts):
    return e0*exp((a*Ts)/(b+Ts))

def gaussian_smoothing(t, x, dx):
    return exp(-((x-t)/dx)**2)

def minus_nabla_F:
    return alpha0 _ alpha1*vbar*dzdx*esat

def integrand():
    return minus_nabla_F *gaussian_smoothing

def expint(x, dx):
    return quad(gaussian_smoothing, x, Inf, args=(x, dx))[0]

vec_expint = vectorize(expint)

