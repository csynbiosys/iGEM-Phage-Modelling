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
S0 = 30.0 #ug/ml(Levin et al., 1977)
D = 0.20 # h-1  #(Levin et al., 1977)
Ki = 6.24e-8 #ml/h (Levin et al., 1977)
b = 98.0 #(Levin et al., 1977)
Km = 4.0 #4 ug/ml(Levin et al., 1977)
Y = 7.40e4 #(Levin et al., 1977)
T = 0.5 #h-1 (Levin et al., 1977)
Xs0 = 1.0e4 # cells/ml starting levels of cells (Levin et al., 1977)
P0 = 0 # particles/ml starting levels of cells (Levin et al., 1977)
q = 0.35 # induction rate (...)
Pt0 = 1.0e6 # particles/ml of temperate phage
Xl = 0 # no lysogenic bacteria present at the start
Xi = 0 # no lytic bacteria present at the start

# Time after start when lytic phage is added
plyt_added = 25.0


dde_camp = p.dde()
dde_camp2 = p.dde()

# Defining the gradient function
# s - state(Xs or P?), c - constant(ddecons), t - time
def ddegrad(s, c, t):
    Xslag = 0.0
    Plag = 0.0

    if (t>c[7]): # if t > T
        Xslag = p.pastvalue(1,t-c[7],0)
        Plag = p.pastvalue(3,t-c[7],0)

    g = array([0.0,0.0,0.0,0.0,0.0,0.0])

    # s[0] = S(t), s[1] = Xs(t), s[2] = Xi(t), s[3] = P(t), s[4] = Xl(t), s[5] = Pt(t)
    # S = D*(S0-S) - Xs*u*S/(Km+S)*(1/Y) - Xi*u*S/(Km+S)*(1/Y) - Xl*u*S/(Km+S)*(1/Y)
    g[0] = c[2]*(c[1]-s[0]) - s[1]*(c[0]*s[0])/(c[5]+s[0])*(1/c[6]) - s[2]*(c[0]*s[0])/(c[5]+s[0])*(1/c[6]) - s[4]*(c[0]*s[0])/(c[5]+s[0])*(1/c[6])
    # Xs = Xs*u*S/(Km+S)*(1/Y) - Ki*Xs*P - Ki*Xs*Pt - D*Xs
    g[1] = s[1]*c[0]*s[0]/(c[5]+s[0]) - c[3]*s[1]*s[3] - c[3]*s[1]*s[4] - c[2]*s[1]
    # Xi = Ki*Xs*P - D*Xi - exp(-D*T)*Ki*Xs(t-T)P(t-T)
    g[2] = 0
    # P = b*exp(-D*T)*Ki*Xs(t-T)P(t-T) - Ki*Xs*P - D*P
    g[3] = 0
    if (t>=plyt_added): # after plyt_added hours lytic phage is added
        # Xi = Ki*Xs*P - D*Xi - exp(-D*T)*Ki*Xs(t-T)P(t-T)
        g[2] = c[3]*s[1]*s[3] - c[2]*s[2] - exp(-c[2]*c[7])*c[3]*Xslag*Plag
        # P = b*exp(-D*T)*Ki*Xs(t-T)P(t-T) - Ki*Xs*P - D*P
        g[3] = c[4]*exp(-c[2]*c[7])*c[3]*Xslag*Plag - c[3]*s[1]*s[3] - c[2]*s[3]
    # Xl = Xl*u*S/(Km+S)*(1/Y) + Ki*Xs*P - q*Xl - -D*Xl
    g[4] = s[4]*(c[0]*s[0])/(c[5]+s[0])*(1/c[6]) + c[3]*s[1]*s[5] -c[8]*s[4]- c[2]*s[4]
    # Pt = bqXl - Ki*Xs*P - D*P
    g[5] = c[4]*c[8]*s[4] - c[3]*s[1]*s[5] - c[2]*s[5]
    return g

# Definte function to store history variables: state and gradient
def ddesthist(g, s, c, t):
    return (s, g)

# Setting constants
ddecons = array([u,S0,D,Ki,b,Km,Y,T,q])
# Setting initial conditions S, Xs, Xi, P, Xl, Pt
S = ddecons[1]
ddeist = array([S, Xs0, Xi, P0, Xl, Pt0])
# Setting a state-scaling array for use in error control when values are very close to 0
ddestsc = array([0,0,0,0,0,0])

# Short version (not used)
#dde_camp.dde(y=ddeist, times=arange(0.0, 250.0, 1.0),
#           func=ddegrad, parms=ddecons,
#           tol=0.000005, dt=1.0, hbsize=10000, nlag=1, ssc=ddestsc)

# Time: 0-plyt_added (Lysogenic phage only)
dde_camp.initproblem(no_vars=6, no_cons=9, nlag=1, nsw=0, t0=0.0, t1=250.0, initstate=ddeist, c=ddecons, otimes= arange(0.0, plyt_added, 0.1), grad=ddegrad, storehistory=ddesthist)
dde_camp.initsolver(tol=0.000005, hbsize=1000, dt=1.0, statescale=ddestsc)
dde_camp.solve()

# Time: plyt_added-250 (Lysogenic + Lytic phage)
# Setting initial conditions = end points from the dde_camp.data
S = dde_camp.data[:, 1][-1]
Xs0 = dde_camp.data[:, 2][-1]
Xi = dde_camp.data[:, 3][-1]
P0 = 1.0e6 # introduction of P lytic
Xl = dde_camp.data[:, 5][-1]
Pt0 = dde_camp.data[:, 6][-1]
# Setting initial conditions S, Xs, Xi, P, Xl, Pt
ddeist = array([S, Xs0, Xi, P0, Xl, Pt0])
# Setting a state-scaling array for use in error control when values are very close to 0
ddestsc = array([0,0,0,0,0,0])

dde_camp2.initproblem(no_vars=6, no_cons=9, nlag=1, nsw=0, t0=0.0, t1=250.0, initstate=ddeist, c=ddecons, otimes= arange(plyt_added, 250.0, 0.1), grad=ddegrad, storehistory=ddesthist)
dde_camp2.initsolver(tol=0.000005, hbsize=1000, dt=1.0, statescale=ddestsc)
dde_camp2.solve()

# Print test values between
print(dde_camp.data[:, 2][-1])
#print(dde_camp2.data[:, 2][0])

plt.plot(concatenate((dde_camp.data[:, 0], dde_camp2.data[:, 0])), concatenate((dde_camp.data[:, 1],dde_camp2.data[:, 1])),  label=r'S')
plt.plot(concatenate((dde_camp.data[:, 0], dde_camp2.data[:, 0])), concatenate((dde_camp.data[:, 2],dde_camp2.data[:, 2])),  label=r'Xs')
plt.plot(concatenate((dde_camp.data[:, 0], dde_camp2.data[:, 0])), concatenate((dde_camp.data[:, 3],dde_camp2.data[:, 3])),  label=r'Xi')
plt.plot(dde_camp2.data[:, 0], dde_camp2.data[:, 4], ":",  label=r'P') # no need to concatenate since previous values = 0
plt.plot(concatenate((dde_camp.data[:, 0], dde_camp2.data[:, 0])), concatenate((dde_camp.data[:, 5],dde_camp2.data[:, 5])),  label=r'X-lys')
plt.plot(concatenate((dde_camp.data[:, 0], dde_camp2.data[:, 0])), concatenate((dde_camp.data[:, 6],dde_camp2.data[:, 6])), ":", label=r'P-temp')

# To test single data sets
# plt.plot(dde_camp.data[:, 0], dde_camp.data[:, 1],  label=r'S')
# plt.plot(dde_camp2.data[:, 0], dde_camp2.data[:, 2],  label=r'Xs')
# plt.plot(dde_camp2.data[:, 0], dde_camp2.data[:, 3],  label=r'Xi')
# plt.plot(dde_camp2.data[:, 0], dde_camp2.data[:, 4],  ":",  label=r'P')
# plt.plot(dde_camp2.data[:, 0], dde_camp2.data[:, 5],  label=r'X-lys')
# plt.plot(dde_camp2.data[:, 0], dde_camp2.data[:, 6],  label=r'P-temp')

plt.legend()
plt.xlabel('Time (hours)')
plt.ylabel('Log concentration (particles/ml)')
plt.yscale('log')
plt.axis([-5,250,1.0e-4,10000000000])
plt.tick_params(
    axis='both', # changes apply to both axis
    labelsize=12) # set new font size
plt.show()
