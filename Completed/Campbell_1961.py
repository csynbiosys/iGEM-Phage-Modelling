import math
from scipy import *
from numpy import *
import PyDDE.pydde as p
import matplotlib.pyplot as plt

# Setting initial values
u = 0.738 # (Levin et al., 1977)
C = 7.0e9 # max carrying capacity (OD600=7)
D = 0.20 # h-1  #(Levin et al., 1977)
Ki = 6.24e-8 #ml/h (Levin et al., 1977)
b = 98.0 # (Levin et al., 1977)
dp = 0.1 # 0... 0.80 h-1 (in full sunlight) (Suttle & Chen, 1992)
T = 0.5 # h-1 (Levin et al., 1977)
Xs0 = 1.0e4 # cells/ml starting leveles of cells (Levin et al., 1977)
P0 = 1.0e6 # particles/ml starting leveles of cells (Levin et al., 1977)

sim_length = 250.0 # set the simulation length time

dde_camp = p.dde()

# Defining the gradient function
# s - state(Xs or P), c - constant(ddecons), t - time
def ddegrad(s, c, t):

    Xslag = 0.0
    Plag = 0.0

    if (t>c[6]): # if t > T
        Xslag = p.pastvalue(0,t-c[6],0)
        Plag = p.pastvalue(1,t-c[6],0)

    g = array([0.0,0.0])

    # s[0] = Xs(t), s[1] = P(t)
    # Xs = u*Xs(t)*(1-Xs(t)/C)-  D*Xs(t)-    Ki*Xs(t)*P(t)
    g[0] = c[0]*s[0]*(1-s[0]/c[1])-c[2]*s[0]-c[3]*s[0]*s[1]
    # P = -D*P(t) -  Ki*Xs(t)*P(t) + b*Ki*Xs(t-T)*P(t-T) - dp*P(t)
    g[1] = -c[2]*s[1] - c[3]*s[0]*s[1] + c[4]*c[3]*Xslag*Plag - c[5]*s[1]

    return g

# Definte function to store history variables: state and gradient
def ddesthist(g, s, c, t):
    return (s, g)

# Setting constants
ddecons = array([u,C,D,Ki,b,dp,T])
# Setting initial conditions Xs and P
ddeist = array([Xs0, P0])
# Setting a state-scaling array for use in error control when values are very close to 0
ddestsc = array([0,0])

# Long version
dde_camp.initproblem(no_vars=2, no_cons=7, nlag=1, nsw=0, t0=0.0, t1=sim_length, initstate=ddeist, c=ddecons, otimes= arange(0.0, sim_length, 0.1), grad=ddegrad, storehistory=ddesthist)
dde_camp.initsolver(tol=0.000005, hbsize=10000, dt=1.0, statescale=ddestsc)
dde_camp.solve()

# Plot figures
plt.style.use('ggplot') # set the global style
p, = plt.plot(dde_camp.data[:, 0], dde_camp.data[:, 2], ":",  label=r'$P$')
xs, = plt.plot(dde_camp.data[:, 0], dde_camp.data[:, 1],  label=r'$X_S$')

f_size = 15 # set font size for plot labels
plt.legend(prop={'size': f_size})
plt.xlabel('Time (hours)', fontsize=f_size)
plt.ylabel('Log concentration (particles/ml)', fontsize=f_size)
plt.yscale('log')
plt.axis([0,sim_length,0,1.0e10])
#plt.text(sim_length*0.8,2.0e9,'$S$= '+str(S0)) # display parameters
plt.tick_params(axis='both', labelsize=f_size)
plt.tight_layout()
plt.show()
plt.savefig('Campbell_PyDDE.pdf')
