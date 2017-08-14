import math
from scipy import *
from numpy import *
import PyDDE.pydde as p
import matplotlib.pyplot as plt

# Setting initial values
u = 0.738 # h-1 (Levin et al., 1977)
S0 = 8.0 # ug/ml(Levin et al., 1977)
D = 0.20 # h-1 (Levin et al., 1977)
Ki = 6.24e-8 #ml/h (Levin et al., 1977)
b = 98.0 # (Levin et al., 1977)
Km = 4.0 # ug/ml (Levin et al., 1977)
Y = 3.85e5 #(Levin et al., 1977)
q = 0.35 # induction rate (Qiu, 2007)
T = 0.0 # no time latency
Xs0 = 1.0e4 # cells/ml starting levels of cells (Levin et al., 1977)
P0 = 1.0e6 # particles/ml starting levels of cells (Levin et al., 1977)

sim_length = 1000.0 # set the simulation length time

dde_camp = p.dde()

# Defining the gradient function
# s - state(Xs,P,Xl,etc), c - constant(ddecons), t - time
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
    # Xs = Xs*u*S/(Km+S) - Ki*Xs*P - D*Xs
    g[1] = s[1]*c[0]*s[0]/(c[5]+s[0]) - c[3]*s[1]*s[3] - c[2]*s[1]
    # Xi = Xl*u*S/(Km+S) + Ki*Xs*P - q*Xl - -D*Xl
    g[2] = s[2]*c[0]*s[0]/(c[5]+s[0]) + c[3]*s[1]*s[3] -c[7]*s[2]- c[2]*s[2]
    # P = bqXl - Ki*Xs*P - D*P
    g[3] = c[4]*c[7]*s[2] - c[3]*s[1]*s[3] - c[2]*s[3]

    return g

# Definte function to store history variables: state and gradient
def ddesthist(g, s, c, t):
    return (s, g)

# Setting constants
ddecons = array([u,S0,D,Ki,b,Km,Y,q,T])
# Setting initial conditions S, Xs, Xi, P
ddeist = array([ddecons[1], Xs0, 0, P0])
# Setting a state-scaling array for use in error control when values are very close to 0
ddestsc = array([0,0,0,0])

# Long version
dde_camp.initproblem(no_vars=4, no_cons=9, nlag=1, nsw=0, t0=0.0, t1=sim_length, initstate=ddeist, c=ddecons, otimes= arange(0.0, sim_length, 0.1), grad=ddegrad, storehistory=ddesthist)
dde_camp.initsolver(tol=0.000005, hbsize=1000, dt=1.0, statescale=ddestsc)
dde_camp.solve()

# Plot figures
plt.style.use('ggplot') # set the global style
pt, = plt.plot(dde_camp.data[:, 0], dde_camp.data[:, 4], "--", label=r'$P_T$')
xs, = plt.plot(dde_camp.data[:, 0], dde_camp.data[:, 2],  label=r'$X_S$')
xl, = plt.plot(dde_camp.data[:, 0], dde_camp.data[:, 3],  label=r'$X_L$')

f_size = 15 # set font size for plot labels
plt.xlabel('Time (hours)', fontsize=f_size)
plt.ylabel('Log concentration (particles/ml)', fontsize=f_size)
plt.yscale('log')
plt.axis([0,sim_length,0.0e-4,1.0e10])
#plt.text(sim_length*0.8,2.0e9,'$S$= '+str(S0)) # display parameters
plt.tick_params(axis='both', labelsize=f_size)

# Plot substrate on the second y axis on top of the preivous figure
# plt2 = plt.twinx()
# plt2.grid(False)
# s, = plt.plot(dde_camp.data[:, 0], dde_camp.data[:, 1], 'black', label=r'$S$')
# plt2.set_ylabel(r'Substrate (${\mu}$g/ml)', fontsize=f_size)
# plt2.set_yticks(linspace(0,S0, 3))
# plt2.tick_params(axis='both', labelsize=f_size)

# Join legends from two separate plots into one
p = [xs,xl,pt]
plt.legend(p, [p_.get_label() for p_ in p],loc='best', fontsize= 'small', prop={'size': f_size})
plt.tight_layout()
#plt.show()
plt.savefig('Qiu8_PyDDE.pdf')
