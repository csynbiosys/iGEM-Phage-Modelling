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

# # Setting initial values (as in Levin 1977)
# u = 7.38 # as in Krysiak-Baltyn, 2016 vs 0.738(Levin et al., 1977)
# S0 = 30.0 #ug/ml(Levin et al., 1977)
# D = 0.20 # h-1  #(Levin et al., 1977)
# Ki = 6.24e-8 #ml/h (Levin et al., 1977)
# b = 98.0 #(Levin et al., 1977)
# Km = 4.0 #4 ug/ml(Levin et al., 1977)
# Y = 7.40e4 #(Levin et al., 1977)
# q = 0.1 # induction rate (...)
# T = 0 # ODL
# Xs0 = 1.0e4 # cells/ml starting levels of cells (Levin et al., 1977)
# P0 = 1.0e6 # particles/ml starting levels of cells (Levin et al., 1977)

# # Setting initial values (as in Levin 1977)
u = 10.0 #0.738 # as in Krysiak-Baltyn, 2016 vs 0.738(Levin et al., 1977)
S0 = 100.0 #ug/ml(Levin et al., 1977)
D = 9.0 #0.20 # h-1  #(Levin et al., 1977)
Ki = 6.24e-8 #ml/h (Levin et al., 1977)
b = 100.0 #(Levin et al., 1977)
Km = 2.0#4.0 #4 ug/ml(Levin et al., 1977)
Y = 0.3#7.40e4 #(Levin et al., 1977)
q = 0.7# 0.35 # induction rate (...)
T = 0 # ODL
Xs0 = 3.789e7 #1.0e4 # cells/ml starting levels of cells (Levin et al., 1977)
P0 = 1.0e6 # particles/ml starting levels of cells (Levin et al., 1977)

dde_camp = p.dde()

# Defining the gradient function
# s - state(Xs or P?), c - constant(ddecons), t - time
def ddegrad(s, c, t):

    Xslag = 0.0
    Plag = 0.0

    if (t>c[7]): # if t > T
        Xslag = p.pastvalue(1,t-c[7],0)
        Plag = p.pastvalue(3,t-c[7],0)

    g = array([0.0,0.0,0.0,0.0])

    # s[0] = S(t), s[1] = Xs(t), s[2] = Xl(t), s[3] = P(t)
    # S = D*(S0-S) - Xs*u*S/(Km+S)*(1/Y) - Xl*u*S/(Km+S)*(1/Y)
    g[0] = c[2]*(c[1]-s[0]) - s[1]*(c[0]*s[0])/(c[5]+s[0])*(1/c[6]) - s[2]*(c[0]*s[0])/(c[5]+s[0])*(1/c[6])
    # Xs = Xs*u*S/(Km+S)*(1/Y) - Ki*Xs*P - D*Xs
    g[1] = s[1]*c[0]*s[0]/(c[5]+s[0]) - c[3]*s[1]*s[3] - c[2]*s[1]
    # Xi = Xl*u*S/(Km+S)*(1/Y) + Ki*Xs*P - q*Xl - -D*Xl
    g[2] = s[2]*(c[0]*s[0])/(c[5]+s[0])*(1/c[6]) + c[3]*s[1]*s[3] -c[7]*s[2]- c[2]*s[2]
    # P = bqXl - Ki*Xs*P - D*P
    g[3] = c[4]*c[7]*s[2] - c[3]*s[1]*s[3] - c[2]*s[3]

    return g

# Definte function to store history variables: state and gradient
def ddesthist(g, s, c, t):
    return (s, g)

# Setting constants
ddecons = array([u,S0,D,Ki,b,Km,Y,q,T])
# Setting initial conditions S, Xs, Xi, P
ddeist = array([ddecons[1], Xs0, 0, P0]) #changed S0 from 100
# Setting a state-scaling array for use in error control when values are very close to 0
ddestsc = array([0,0,0,0])

# Short version (not used)
#dde_camp.dde(y=ddeist, times=arange(0.0, 250.0, 1.0),
#           func=ddegrad, parms=ddecons,
#           tol=0.000005, dt=1.0, hbsize=10000, nlag=1, ssc=ddestsc)

# Long version
dde_camp.initproblem(no_vars=4, no_cons=9, nlag=1, nsw=0, t0=0.0, t1=250.0, initstate=ddeist, c=ddecons, otimes= arange(0.0, 250.0, 0.1), grad=ddegrad, storehistory=ddesthist)

dde_camp.initsolver(tol=0.000005, hbsize=1000, dt=1.0, statescale=ddestsc)

dde_camp.solve()

print(dde_camp.data)

#plt.plot(dde_camp.data[:, 0], dde_camp.data[:, 1],  label=r'S')
plt.plot(dde_camp.data[:, 0], dde_camp.data[:, 2],  label=r'Xs')
#plt.plot(dde_camp.data[:, 0], dde_camp.data[:, 3],  label=r'Xl')
plt.plot(dde_camp.data[:, 0], dde_camp.data[:, 4],  label=r'P')
#plt.plot(dde_camp.data[:, 0], dde_camp.data[:, 1],  label=r'S')
plt.legend()
plt.xlabel('Time (hours)')
plt.ylabel('Log concentration (particles/ml)')
plt.yscale('log')
plt.axis([-5,250,-100,10000000000])
plt.tick_params(
    axis='both', # changes apply to both axis
    labelsize=12) # set new font size
plt.show()
