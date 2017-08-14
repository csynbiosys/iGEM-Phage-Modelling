import math
from scipy import *
from numpy import *
import PyDDE.pydde as p
import matplotlib.pyplot as plt

# Setting initial values
u = 0.738 # h-1 (Levin et al., 1977)
S0 = 3.0e1 # ug/ml(Levin et al., 1977)
Km = 4.0 # ug/ml (Levin et al., 1977)
Y = 3.85e5 #(Levin et al., 1977)
Xs0 = 2.3e9 # cells/ml starting levels of cells (Levin et al., 1977)
a = 0.22 # death related constant (Kumada et al., 1985)
k = -2.0e-4 # death related constant (Kumada et al., 1985)
C = 7.0e9 # max carrying capacity (OD600=7)

sim_length = 800.0 # set the simulation length time

dde_camp = p.dde()

# Defining the gradient function
# s - state(Xs,P,Xl,etc), c - constant(ddecons), t - time
def ddegrad(s, c, t):

    g = array([0.0,0.0])

    # s[0] = S(t), s[1] = Xs(t)
    # S = Xs*u*S/(Km+S)*(1/Y)
    g[0] = - s[1]*(c[0]*s[0])/(c[5]+s[0])*(1/c[6])
    # X = k*X**(a+1)
    g[1] = (c[2]*s[1]**c[3])*s[1]
    return g

# Definte function to store history variables: state and gradient
def ddesthist(g, s, c, t):
    return (s, g)

# Setting constants
ddecons = array([u,S0,k,a,Xs0,Km,Y,C])
# Setting initial conditions S, Xs, Xi, P
ddeist = array([ddecons[1], Xs0]) #changed S0 from 100
# Setting a state-scaling array for use in error control when values are very close to 0
ddestsc = array([0,0])

# Long version
dde_camp.initproblem(no_vars=2, no_cons=8, nlag=1, nsw=0, t0=0.0, t1=sim_length, initstate=ddeist, c=ddecons, otimes= arange(0.0, sim_length, 20.0), grad=ddegrad, storehistory=ddesthist)

dde_camp.initsolver(tol=0.000005, hbsize=1000, dt=1.0, statescale=ddestsc)

dde_camp.solve()

print(dde_camp.data)

f_size = 14 # set size for plot labels
#plt.plot(dde_camp.data[:, 0], dde_camp.data[:, 1],  label=r'$S$')
plt.plot(dde_camp.data[:, 0], dde_camp.data[:, 2],  label=r'$X_S$')

plt.legend(prop={'size': f_size})
plt.xlabel('Time (hours)', fontsize=f_size)
plt.ylabel('Log concentration (particles/ml)', fontsize=f_size)
plt.yscale('log')
plt.axis([-5,sim_length, 1.0e0,1.0e10])
plt.tick_params(
    axis='both', # changes apply to both axis
    labelsize=f_size) # set new font size
plt.tight_layout()
plt.show()
