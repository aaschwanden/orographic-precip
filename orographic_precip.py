#!/usr/bin/env python

import numpy as np
from scipy.integrate import odeint, quad
import pylab as plt



def f(x):
    "surface elevation"
    return 1./(np.sqrt(2*np.pi)) * np.exp(-0.5*(x - 5.0)**2)

def g(x):
    a = 4
    return f(x) * quad(a * f(x), -np.inf, np.inf)[0]

def z(x, a=1, z0=0):
    '''
    Right-skewed topography
    '''
    def f(x):
        return 1./(np.sqrt(2*np.pi)) * np.exp(-0.5*(x - 5.0)**2)
    def g(x):
        return a * 1./(np.sqrt(2*np.pi)) * np.exp(-0.5*(x - 5.0)**2)

    return f(x) / z0 * quad(g, -np.inf, np.inf)



# def z(x):
#     "surface elevation"
#     return 0.5 * (3.0 * np.exp(-(x - 5.0)**2) + 2.0 * np.exp(-(x - 9.0)**2) + 1.0 * np.exp(-(x - 13.0)**2))

def h(x):
    "ice elevation"
    return 0.25 * (3.0 * np.exp(-(x - 4.0)**2) + 2.0 * np.exp(-(x - 8.0)**2) + 1.0 * np.exp(-(x - 12.0)**2))

def es(T):
    "saturation vapor pressure"
    return 6.1094 * np.exp((17.625 * T) / (T + 243.04))

def Tm(x):
    "air temperature"
    return 10.0 - 6.0 * z(x)    # lapse rate

def Ti(x):
    "air temperature"
    return 10.0 - 6.0 * (z(x) + h(x))    # lapse rate

def windwardm(x):
    dx = 0.01
    dz = z(x + dx) - z(x)
    if dz / dx > 0:
        return 1
    else:
        return 0

def windwardi(x):
    dx = 0.01
    dz = z(x + dx) + h(x + dx) - (z(x) + h(x))
    if dz / dx > 0:
        return 1
    else:
        return 0
    
alpha0 = 0.01
alpha1 = 0.1

def pm(y, x):
    "right hand side"
    dx = 0.01
    dz = z(x + dx) - z(x)
    return -y * es(Tm(x)) * (alpha0 + alpha1 * windwardm(x) * v * (dz/dx))

def pi(y, x):
    "right hand side"
    dx = 0.01
    dz = z(x + dx) + h(x + dx) - (z(x) + h(x))
    return -y * es(Ti(x)) * (alpha0 + alpha1 * windwardi(x) * v * (dz/dx))


x = np.linspace(0, 20, 1001)
y = z(x)
yi = z(x) + h(x) 


fig = plt.figure()
v = 0.5
Pm = np.squeeze(odeint(pm, 1.0, x))
dPm = - np.diff(Pm, axis=0) / (x[1] - x[0])
Pi = np.squeeze(odeint(pi, 1.0, x))
dPi = - np.diff(Pi, axis=0) / (x[1] - x[0])

axU = fig.add_subplot(211)
axU.plot(x[1:], dPm, color='blue')
axU.plot(x[1:], dPi, color='red')
axL = fig.add_subplot(212)
axL.plot(x[1:], y[1:], color='blue')
axL.plot(x[1:], yi[1:], color='red')


# for k in range(1,6):
#     v = 0.5 * k
#     P = np.squeeze(odeint(p, 1.0, x))
#     dP = - np.diff(P, axis=0) / (x[1] - x[0])

#     plt.subplot(5, 1, k)
#     plt.hold(True)
#     plt.plot(x[1:], dP, lw=1, color="blue", label="precipitation for v=%f" % v)
#     plt.plot(x[1:], y[1:] + dP, "--", lw=1, color="red", label="precipitation over geometry")
#     plt.plot(x, y, lw=2, color="black", label="top surface elevation")
#     plt.legend(loc='best')

plt.grid(True)
plt.show()
