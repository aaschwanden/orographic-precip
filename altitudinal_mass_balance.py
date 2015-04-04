#!/usr/bin/env python
# Copyright (C) 2015 Andy Aschwanden

# nosetests altitudinal_mass_balance.py

import numpy as np
import pylab as plt

# Set up constants in SI units
L = 60e3  # m
nx =3e2  #
dx = 30e3  # m

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

    
def M_s(T_0, z, lat):

    gamma_1 = 1.0

    return -gamma_1 * (T_s(T_0, z, lat) - 0.5)

def T_s(T_0, z, lat):

    alpha = 0.01   # degC / m

    return T_0 - alpha * z - 0.7576 * lat

z = gaussian_bump(x, 2e3, 0, L/2)
z = np.linspace(0, 3000, 301)


T_0 = 46  # There must be a typo somewhere in the manuscript, the latitude
          # dependence is way too strong.

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(M_s(T_0., z, 43.5), z)

plt.show()
