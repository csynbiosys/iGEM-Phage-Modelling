#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 12 14:10:03 2017

@author: antonpuzorjov
"""

import math
from scipy import *
import numpy as np
#import PyDDE.pydde as p
import matplotlib.pyplot as plt
from jitcdde import provide_basic_symbols, jitcdde

#DDE example from Solv95 distribution.

#This model is a model for Nicholson's (1954) blowflies, as given by Gurney and Nisbet (1981)
P = 10.0
A0 = 300.0
d = 0.25
T = 12.0

# A() = y(0); A(t-T) = y(0, t-T)
t, y = provide_basic_symbols()

f = [P*y(0,t-T)*exp(-y(0,t-T)/A0)-d*y(0)] # ERROR AttributeError: 'Mul' object has no attribute 'exp'

# initialising the integrator
DDE = jitcdde(f)

# enter initial conditions
#N0 = 0.1
#No0 = 10
DDE.add_past_point(-1.0, [1.0], [0.0])
DDE.add_past_point( 0.0, [1.0], [0.0])

# short pre-integration to take care of discontinuities
DDE.step_on_discontinuities()

# create timescale
stoptime = 1000.0
numpoints = 10000
times = DDE.t + linspace(0, stoptime, numpoints)

# integrating
data = []
for time in times:
    data.append(DDE.integrate(time))
    
print(data)



#
#
#dde_eg = p.dde()
#
#def ddegrad(s, c, t):
#    alag = 0.0
#    if (t>c[0]):
#        alag = p.pastvalue(0,t-c[0],0)
#    return array( [ c[2]*alag*exp(-alag/c[3])-c[1]*s[0] ] )
#
#def ddesthist(g, s, c, t):
#    return (s, g)
#
#ddecons = array([12.0,0.25,10.0,300.0,100.0])
#ddeist = array([ddecons[4]])
#ddestsc = array([0])
#
## Short version
#dde_eg.dde(y=ddeist, times=arange(0.0, 300.0, 1.0),
#           func=ddegrad, parms=ddecons,
#           tol=0.000005, dt=1.0, hbsize=1000, nlag=1, ssc=ddestsc)
#
##dde_eg.initproblem(no_vars=1, no_cons=5, nhv=1, nlag=1, nsw=0, no_otimes=301, t0=0.0, t1=300.0, initstate=ddeist, c=ddecons, otimes= arange(0.0, 300.0, 1.0), grad=ddegrad, storehistory=ddesthist)
#
#
## print(dde_eg.data)
#
#plt.plot(dde_eg.data[:, 0], dde_eg.data[:, 1],  label=r'blowflies')
#plt.legend()
#plt.xlabel('time t')
#plt.ylabel('population size')
##plt.yscale('log')
#plt.axis([0,300,10,4000])
#plt.tick_params(
#    axis='both', # changes apply to both axis
#    labelsize=12) # set new font size
