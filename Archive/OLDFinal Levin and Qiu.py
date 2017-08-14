import math
from scipy import *
from numpy import *
import PyDDE.pydde as p
import matplotlib.pyplot as plt

sim_length = 25.0 # set the simulation length time
plyt_added = 5.0 # time after start when lytic phage is added

# Setting initial values
u = 0.738 # h-1 (Levin et al., 1977)
S0 = 100.0 # ug/ml(Levin et al., 1977)
D = 0.20 # h-1  HERE USED AS COMMON DEATH RATE FOR BACTERIA
Ki = 6.24e-8 # ml/h (Levin et al., 1977)
b = 98.0 # (Levin et al., 1977)
Km = 4.0 # ug/ml(Levin et al., 1977)
Y = 7.40e4 # (Levin et al., 1977)
T = 0.5 #h-1 (Levin et al., 1977)
Xs0 = 1.0e6 # cells/ml starting levels of cells (Levin et al., 1977)
P0 = 0 # particles/ml starting levels of cells (Levin et al., 1977)
q = 0.001 # induction rate (...)
Pt0 = 1.0e6 # particles/ml of temperate phage
Xl = 0 # no lysogenic bacteria present at the start
Xi = 0 # no lytic bacteria present at the start


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
    # S = - Xs*u*S/(Km+S)*(1/Y) - Xi*u*S/(Km+S)*(1/Y) - Xl*u*S/(Km+S)*(1/Y)
    g[0] = - s[1]*(c[0]*s[0])/(c[5]+s[0])*(1/c[6]) - s[2]*(c[0]*s[0])/(c[5]+s[0])*(1/c[6]) - s[4]*(c[0]*s[0])/(c[5]+s[0])*(1/c[6])
    # Xs = Xs*u*S/(Km+S)*(1/Y) - Ki*Xs*P - Ki*Xs*Pt - D*Xs
    g[1] = s[1]*c[0]*s[0]/(c[5]+s[0]) - c[3]*s[1]*s[3] - c[3]*s[1]*s[5] - c[2]*s[1]
    # Xi = Ki*Xs*P - D*Xi - Ki*Xs(t-T)P(t-T)
    g[2] = 0
    # P = b*exp(-D*T)*Ki*Xs(t-T)P(t-T) - Ki*Xs*P
    g[3] = 0
    if (t>plyt_added): # after plyt_added hours lytic phage is added
        # Xi = Ki*Xs*P - D*Xi - Ki*Xs(t-T)P(t-T)
        g[2] = c[3]*s[1]*s[3] - c[2]*s[2] - c[3]*Xslag*Plag
        # P = b*Ki*Xs(t-T)P(t-T) - Ki*Xs*P - Ki*Xl*P
        g[3] = c[4]*c[3]*Xslag*Plag - c[3]*s[1]*s[3] - c[3]*s[4]*s[3]
    # Xl = Xl*u*S/(Km+S)*(1/Y) + Ki*Xs*Pt - q*Xl - -D*Xl   # q*Xl induction rate still exists but set for very low value
    g[4] = s[4]*(c[0]*s[0])/(c[5]+s[0])*(1/c[6]) + c[3]*s[1]*s[5] -c[8]*s[4]- c[2]*s[4]
    # Pt = bqXl - Ki*Xs*Pt
    g[5] = c[4]*c[8]*s[4] - c[3]*s[1]*s[5]
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
dde_camp.initproblem(no_vars=6, no_cons=9, nlag=1, nsw=0, t0=0.0, t1=plyt_added, initstate=ddeist, c=ddecons, otimes= arange(0.0, plyt_added, 0.1), grad=ddegrad, storehistory=ddesthist)
dde_camp.initsolver(tol=0.000005, hbsize=10000, dt=1.0, statescale=ddestsc)
dde_camp.solve()

# Time: plyt_added-250 (Lysogenic + Lytic phage)
# Setting initial conditions
S = dde_camp.data[:, 1][-1]
Xs0 = dde_camp.data[:, 2][-1]
Xi = dde_camp.data[:, 3][-1]
P0 = 1.0e4 # introduction of P lytic
Xl = dde_camp.data[:, 5][-1]
Pt0 = dde_camp.data[:, 6][-1]
# Setting initial conditions S, Xs, Xi, P, Xl, Pt
ddeist = array([S, Xs0, Xi, P0, Xl, Pt0])
# Setting a state-scaling array for use in error control when values are very close to 0
ddestsc = array([0,0,0,0,0,0])

dde_camp2.initproblem(no_vars=6, no_cons=9, nlag=1, nsw=0, t0=plyt_added, t1=sim_length, initstate=ddeist, c=ddecons, otimes= arange(plyt_added, sim_length, 0.1), grad=ddegrad, storehistory=ddesthist)
dde_camp2.initsolver(tol=0.000005, hbsize=10000, dt=1.0, statescale=ddestsc)
dde_camp2.solve()

# Plot figures
plt.style.use('ggplot') # set the global style
xs, = plt.plot(concatenate((dde_camp.data[:, 0], dde_camp2.data[:, 0])), concatenate((dde_camp.data[:, 2],dde_camp2.data[:, 2])),  label=r'$X_S$')
xi, = plt.plot(concatenate((dde_camp.data[:, 0], dde_camp2.data[:, 0])), concatenate((dde_camp.data[:, 3],dde_camp2.data[:, 3])),  label=r'$X_I$')
p, = plt.plot(dde_camp2.data[:, 0], dde_camp2.data[:, 4], "r:",  label=r'$P$') # no need to concatenate since previous values = 0
xl, = plt.plot(concatenate((dde_camp.data[:, 0], dde_camp2.data[:, 0])), concatenate((dde_camp.data[:, 5],dde_camp2.data[:, 5])),  label=r'$X_L$')
pt, = plt.plot(concatenate((dde_camp.data[:, 0], dde_camp2.data[:, 0])), concatenate((dde_camp.data[:, 6],dde_camp2.data[:, 6])), "r--", label=r'$P_T$')

f_size = 15 # set font size for plot labels
plt.xlabel('Time (hours)', fontsize=f_size)
plt.ylabel('Log concentration (particles/ml)', fontsize=f_size)
plt.yscale('log')
plt.axis([0,sim_length,1.0e-4,1.0e10])
#plt.text(sim_length*0.8,2.0e9,'$S$= '+str(S0)) # display parameters
plt.tick_params(axis='both', labelsize=f_size)

# Plot substrate on the second y axis on top of the preivous figure
plt2 = plt.twinx()
plt2.grid(False)
s, = plt.plot(concatenate((dde_camp.data[:, 0], dde_camp2.data[:, 0])), concatenate((dde_camp.data[:, 1],dde_camp2.data[:, 1])), 'black',  label=r'$S$')
plt2.set_ylabel(r'Substrate (${\mu}$g/ml)', fontsize=f_size)
plt2.set_yticks(linspace(0,S0, 3))
plt2.tick_params(axis='both', labelsize=f_size)

# Join legends from two separate plots into one
p = [xs,xi,p,xl,pt,s]
plt.legend(p, [p_.get_label() for p_ in p],loc='best', fontsize= 'small', prop={'size': f_size})
plt.tight_layout()
plt.show()
#plt.savefig('Model Batch Final.pdf')
