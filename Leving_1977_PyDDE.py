#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 13 19:23:17 2017

@author: antonpuzorjov
"""

import math
#%matplotlib inline
from scipy import *
from numpy import *
import PyDDE.pydde as p
import matplotlib.pyplot as plt


# Standard parameters
# u = 7.38 (Levin et al., 1977)
# Ki = 10^(-7) ml/h (Levin et al., 1977)
# b = 90 (Levin et al., 1977)
# D = 0.20 h-1  (Levin et al., 1977)
# dp = 0 ... 0.80 h-1 (in full sunlight) (Suttle & Chen, 1992)
# C = 3.5*10^9 (max carrying capacity (OD600=7))
# T = 0.5 h-1 (Levin et al., 1977)
# Xs0 = 2.25*10^4 # cells/ml.
# P0 = 5*(10^6) # particles/ml starting levels of cells (Levin et al., 1977)

# Setting initial values
u = 7.38 #(Levin et al., 1977)
S0 = 30 #ug/ml(Levin et al., 1977)
D = 0.20 # h-1  #(Levin et al., 1977)
Ki = 6.24*10**(-8) #ml/h (Levin et al., 1977)
b = 98 #(Levin et al., 1977)
Km = 4 #4 ug/ml(Levin et al., 1977)
Y = 7.40*(10**4) #(Levin et al., 1977)
T = 0.5 #h-1 (Levin et al., 1977)
Xs0 = 2.25*(10**4) # cells/ml starting levels of cells (Levin et al., 1977)
P0 = 5*(10**6) # particles/ml starting levels of cells (Levin et al., 1977)

dde_camp = p.dde()

# Defining the gradient function
# s - state(Xs or P?), c - constant(ddecons), t - time
def ddegrad(s, c, t):
    g = array([0.0,0.0,0.0,0.0])
    
    # s[0] = S(t), s[1] = Xs(t), s[2] = Xi(t), s[3] = P(t)
    # S = D*(S0-S) - Xs*u*S/(Km+S)*(1/Y) - Xi*u*S/(Km+S)*(1/Y)
    g[0] = c[2]*(c[1]-s[0]) - s[1]*(c[0]*s[0])/(c[5]+s[0])*(1/c[6]) - s[2]*(c[0]*s[0])/(c[5]+s[0])*(1/c[6])
    # Xs = Xs*u*S/(Km+S)*(1/Y) - Ki*Xs*P - D*Xs
    g[1] = s[1]*c[0]*s[0]/(c[5]+s[0]) - c[3]*s[1]*s[3] - c[2]*s[1]
    # Xi = Ki*Xs*P - D*Xi - exp(-D*T)*Ki*Xs(t-T)P(t-T) # without delay
    g[2] = c[5]*s[1]*s[3] - c[2]*s[2] - exp(-c[2]*c[7])*c[3]*s[1]*s[3]
    # P = b*exp(-D*T)*Ki*Xs(t-T)P(t-T) - Ki*Xs*P - D*P  # without delay
    g[3] = c[4]*exp(-c[2]*c[7])*c[3]*s[1]*s[3] - c[3]*s[1]*s[3] - c[2]*s[3]
    
    if (t>c[7]): # if t > T
        # Xi with delay
        g[2] = c[5]*s[1]*s[3] - c[2]*s[2] - exp(-c[2]*c[7])*c[3]*p.pastvalue(1,t-c[7],0)*p.pastvalue(3,t-c[7],0)
        # P with delay
        g[3] = c[4]*exp(-c[2]*c[7])*c[3]*p.pastvalue(1,t-c[7],0)*p.pastvalue(3,t-c[7],0) - c[3]*s[1]*s[3] - c[2]*s[3]
    return g

# Definte function to store history variables: state and gradient
def ddesthist(g, s, c, t):
    return (s, g)

# Setting constants
ddecons = array([u,S0,D,Ki,b,Km,Y,T,Xs0,P0])
# Setting initial conditions S, Xs, Xi, P
ddeist = array([100, ddecons[8], 0, ddecons[9]])
# Setting a state-scaling array for use in error control when values are very close to 0
ddestsc = 0*ddeist

# Short version
dde_camp.dde(y=ddeist, times=arange(0.0, 250.0, 1.0),
           func=ddegrad, parms=ddecons,
           tol=0.000005, dt=1.0, hbsize=10000, nlag=1, ssc=ddestsc)

# Long version (failed to work)
#dde_eg.initproblem(no_vars=1, no_cons=5, nhv=1, nlag=1, nsw=0, no_otimes=301, t0=0.0, t1=300.0, initstate=ddeist, c=ddecons, otimes= arange(0.0, 300.0, 1.0), grad=ddegrad, storehistory=ddesthist)


print(dde_camp.data)

#plt.plot(dde_camp.data[:, 0], dde_camp.data[:, 1],  label=r'S')
plt.plot(dde_camp.data[:, 0], dde_camp.data[:, 2],  label=r'Xs')
plt.plot(dde_camp.data[:, 0], dde_camp.data[:, 3],  label=r'Xi')
plt.plot(dde_camp.data[:, 0], dde_camp.data[:, 4],  label=r'P')
plt.legend()
plt.xlabel('time t')
plt.ylabel('population size')
plt.yscale('log')
plt.axis([-5,25,-1000000,10000000000])
plt.tick_params(
    axis='both', # changes apply to both axis
    labelsize=12) # set new font size