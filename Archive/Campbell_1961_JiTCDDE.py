#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 13 12:48:20 2017

@author: antonpuzorjov
"""

import math
#%matplotlib inline
from scipy import *
from numpy import *
import PyDDE.pydde as p
import matplotlib.pyplot as plt
from jitcdde import provide_basic_symbols, jitcdde

# Standard parameters
# u = 7.38 (Levin et al., 1977)
# Ki = 10^(-7) ml/h (Levin et al., 1977)
# b = 90 (Levin et al., 1977)
# D = 0.20 h-1  (Levin et al., 1977)
# dp = 0 ... 0.80 h-1 (in full sunlight) (Suttle & Chen, 1992)
# C = 3.5*10^9 (max carrying capacity (OD600=7))
# T = 0.5 h-1 (Levin et al., 1977)
# Xs0 = 2.25*10^4 # cells/ml.
# P0 = 5*(10^6) # particles/ml starting leveles of cells (Levin et al., 1977)

# Setting initial values
u = 7.38 #(Levin et al., 1977)
C = 3.5*(10**9) #(max carrying capacity (OD600=7))
D = 0.20 # h-1  #(Levin et al., 1977)
Ki = 10**(-7) #ml/h (Levin et al., 1977)
b = 90 #(Levin et al., 1977)
dp = 0.5#... 0.80 h-1 (in full sunlight) (Suttle & Chen, 1992)
T = 0.5 #h-1 (Levin et al., 1977)
Xs0 = 2.25*(10**4) # cells/ml starting leveles of cells (Levin et al., 1977)
P0 = 5*(10**6) # particles/ml starting leveles of cells (Levin et al., 1977)

# Xs() = y(0); Xs(t-T) = y(0, t-T)
# P() = y(1); X(t-T) = y(1, t-T)
t, y = provide_basic_symbols()

f = [u*y(0)*(1-y(0)/C)-D*y(0)-Ki*y(0)*y(1),
     -D*y(1)-Ki*y(0)*y(1)+b*Ki*y(0, t-T)*y(1, t-T)-dp*y(1)
     ]

# initialising the integrator
DDE = jitcdde(f)

# enter initial conditions
#N0 = 0.1
#No0 = 10
DDE.add_past_point(-1.0, [Xs0,P0], [0.0,0.0])
DDE.add_past_point(0.0, [Xs0,P0], [0.0,0.0])

DDE.set_integration_parameters(atol=10**(-12)) # changing integration parameters

# short pre-integration to take care of discontinuities
DDE.step_on_discontinuities()

# create timescale
stoptime =40.0 # Does not integrate past ~30
numpoints = 100
times = DDE.t + linspace(0, stoptime, numpoints)
size = times.size

# integrating
#data = []#empty((size,3))
data=[]
for time in times:
    data = append(data, [[time, DDE.integrate(time)[0], DDE.integrate(time)[1]]])

data = data.reshape(size,3)    
print(data)



plt.plot(data[:,0], data[:,1],  label=r'Xs')
plt.plot(times, data[:, 2],  label=r'P')
plt.legend()
plt.xlabel('time t')
plt.ylabel('population size')
plt.yscale('log')
plt.axis([0,50,0,10000000000])
plt.tick_params(
    axis='both', # changes apply to both axis
    labelsize=12) # set new font size

#
#
#dde_camp = p.dde()
#
## Defining the gradient function
## s - state(Xs or P?), c - constant(ddecons), t - time?
#def ddegrad(s, c, t):
#    g = array([0.0,0.0])
#    # u*Xs(t)*(1-Xs(t)/C)-  D*Xs(t)-    Ki*Xs(t)*P(t)
#    # s[0] = Xs(t), s[1] = P(t)
#    g[0] = c[0]*s[0]*(1-s[0]/c[1])-c[2]*s[0]-c[3]*s[0]*s[1] # dXs/dt
#    # -D*P(t)-  Ki*Xs(t)*P(t)+  b*Ki*Xs(t-T)*P(t-T)-    dp*P(t)
#    g[1] = -c[2]*s[1]-c[3]*s[0]*s[1]+c[4]*c[3]*s[0]*s[1]-c[5]*s[1] # dP/dt without delay
#    if (t>c[6]): # if t > T
#        g[1] = -c[2]*s[1]-c[3]*s[0]*s[1]+c[4]*c[3]*p.pastvalue(0,t-c[6],0)*p.pastvalue(1,t-c[6],0)-c[5]*s[1] # dP/dt without delay
#    return g
#
## Definte function to store history variables: state and gradient
#def ddesthist(g, s, c, t):
#    return (s, g)
#
## Setting constants
#ddecons = array([u,C,D,Ki,b,dp,T,Xs0,P0])
## Setting initial conditions
#ddeist = array([ddecons[7], ddecons[8]])
## Setting a state-scaling array for use in error control when values are very close to 0
#ddestsc = 0*ddeist
#
## Short version
#dde_camp.dde(y=ddeist, times=arange(0.0, 250.0, 1.0),
#           func=ddegrad, parms=ddecons,
#           tol=0.000005, dt=1.0, hbsize=1000, nlag=1, ssc=ddestsc)
#
## Long version (failed to work)
##dde_eg.initproblem(no_vars=1, no_cons=5, nhv=1, nlag=1, nsw=0, no_otimes=301, t0=0.0, t1=300.0, initstate=ddeist, c=ddecons, otimes= arange(0.0, 300.0, 1.0), grad=ddegrad, storehistory=ddesthist)
#
#
#print(dde_camp.data)
#
#plt.plot(dde_camp.data[:, 0], dde_camp.data[:, 1],  label=r'Xs')
#plt.plot(dde_camp.data[:, 0], dde_camp.data[:, 2],  label=r'P')
#plt.legend()
#plt.xlabel('time t')
#plt.ylabel('population size')
#plt.yscale('log')
#plt.axis([-5,250,-1000000,10000000000])
#plt.tick_params(
#    axis='both', # changes apply to both axis
#    labelsize=12) # set new font size