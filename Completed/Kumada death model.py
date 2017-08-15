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
Xs0 = 1.0e4 # cells/ml starting levels of cells (Levin et al., 1977)
a = 0.22 # death related constant (Kumada et al., 1985)
k = -2.0e-4 # death related constant (Kumada et al., 1985)
C = 7.0e9 # max carrying capacity (OD600=7)

sim_length = 2500.0 # set the simulation length time

dde_camp = p.dde()

# Defining the gradient function
# s - state(Xs,P,Xl,etc), c - constant(ddecons), t - time
def ddegrad(s, c, t):

    g = array([0.0,0.0])

    # s[0] = S(t), s[1] = Xs(t)
    # S = Xs*u*S/(Km+S)*(1/Y)
    g[0] = - s[1]*(c[0]*s[0])/(c[5]+s[0])*(1/c[6])
    # check if S is low or cells are over C
    if (s[0] < 1/c[6] or s[1]>c[7]):
        # cells start to die
        # X = k*X**(a+1)
        g[1] = c[2]*s[1]**(c[3]+1)
    else:
        # cells grow normally
        # X = Xs*u*S/(Km+S)
        g[1] = s[1]*c[0]*s[0]/(c[5]+s[0])

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
dde_camp.initproblem(no_vars=2, no_cons=8, nlag=1, nsw=0, t0=0.0, t1=sim_length, initstate=ddeist, c=ddecons, otimes= arange(0.0, sim_length, 0.1), grad=ddegrad, storehistory=ddesthist)
dde_camp.initsolver(tol=0.000005, hbsize=1000, dt=1.0, statescale=ddestsc)
dde_camp.solve()


# Plot figures
plt.style.use('ggplot') # set the global style
xs, = plt.plot(dde_camp.data[:, 0], dde_camp.data[:, 2], '#398aba',  label=r'$X_S$')

f_size = 15 # set font size for plot labels
plt.xlabel('Time (hours)', fontsize=f_size)
plt.ylabel('Log concentration (particles/ml)', fontsize=f_size)
plt.yscale('log')
plt.axis([0,sim_length,0,1.0e10])
plt.tick_params(axis='both', labelsize=f_size)

# Plot substrate on the second y axis on top of the preivous figure
plt2 = plt.twinx()
plt2.grid(False)
s, = plt.plot(dde_camp.data[:, 0], dde_camp.data[:, 1],  label=r'$S$')
plt2.set_ylabel(r'Substrate (${\mu}$g/ml)', fontsize=f_size)
plt2.set_yticks(linspace(0,S0, 3))
plt2.tick_params(axis='both', labelsize=f_size)

# Join legends from two separate plots into one
p = [xs,s]
plt.legend(p, [p_.get_label() for p_ in p],loc='best', fontsize= 'small', prop={'size': f_size})
plt.tight_layout()
plt.show()
#plt.savefig('Kumada_with_growth_1985.pdf')
