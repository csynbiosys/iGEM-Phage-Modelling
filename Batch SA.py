import math
from scipy import *
from numpy import *
import PyDDE.pydde as p
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick

# Nominal values
u_nominal = 0.738 # h-1 (Levin et al., 1977)
S0_nominal = 30.0 # ug/ml(Levin et al., 1977)
D = 0.0 # h-1 not used here
Ki_nominal = 6.24e-8 # ml/h (Levin et al., 1977)
b_nominal = 98.0 # (Levin et al., 1977)
Km_nominal = 4.0 # ug/ml(Levin et al., 1977)
Y_nominal = 3.85e5 # (Levin et al., 1977)
T_nominal = 0.5 #h-1 (Levin et al., 1977)
Xs0_nominal = 1.0e4 # cells/ml starting levels of cells (Levin et al., 1977)
P0_nominal = 0 # particles/ml starting levels of cells (Levin et al., 1977)
q_nominal = 0.35 # induction rate (Qiu, 2007)
Pt0_nominal = 1.0e6 # particles/ml of temperate phage
# Extended parameters for Xi, Xl and Pt
ui_nominal = 0.738 # for lytically infected bacteria
ul_nominal = 0.738 # for lysogenic bacteria
Kmi_nominal = 4.0 # for lytically infected bacteria
Kml_nominal = 4.0 # for lysogenic bacteria
Yi_nominal = 3.85e5 # for lytically infected bacteria
Yl_nominal = 3.85e5 # for lysogenic bacteria
Tt_nominal = 0.5 # for temperate phage
Kit_nominal = 6.24e-8 # for temperate phage
bt_nominal = 98.0 # for temperate phage

sim_length_nominal = 20.0 # set the simulation length time
plyt_added_nominal = 5.0 # time after start when lytic phage is added

# Setting initial values
u = u_nominal # h-1 (Levin et al., 1977)
S0 = S0_nominal # ug/ml(Levin et al., 1977)
D = 0.0 # h-1 not used here
Ki = Ki_nominal # ml/h (Levin et al., 1977)
b = b_nominal # (Levin et al., 1977)
Km = Km_nominal # ug/ml(Levin et al., 1977)
Y = Y_nominal # (Levin et al., 1977)
T = T_nominal #h-1 (Levin et al., 1977)
Xs0 = Xs0_nominal # cells/ml starting levels of cells (Levin et al., 1977)
P0 = P0_nominal # particles/ml starting levels of cells (Levin et al., 1977)
q = q_nominal # induction rate (Qiu, 2007)
Pt0 = Pt0_nominal # particles/ml of temperate phage
Xl = 0 # no lysogenic bacteria present at the start
Xi = 0 # no lytic bacteria present at the start

# Extended parameters for Xi, Xl and Pt
ui = ui_nominal # for lytically infected bacteria
ul = ul_nominal # for lysogenic bacteria
Kmi = Kmi_nominal # for lytically infected bacteria
Kml = Kml_nominal # for lysogenic bacteria
Yi = Yi_nominal # for lytically infected bacteria
Yl = Yl_nominal # for lysogenic bacteria
Tt = Tt_nominal  # for temperate phage
Kit = Kit_nominal # for temperate phage
bt = bt_nominal # for temperate phage

sim_length = sim_length_nominal # set the simulation length time
plyt_added = plyt_added_nominal # time after start when lytic phage is added



def dde_sa (parameter, min, max, step):
    result = {} # dictionary to save {percentage : [Xs_extinct, xixs_ratio]}
    for percentage in linspace(min,max,1):
        if parameter == 'u':
            global u
            u = u_nominal*percentage/100
            print('Percentage is ' + str(percentage))
            print('New u= ' + str(u))

        # Defining the gradient function
        # s - state(Xs or P?), c - constant(ddecons), t - time
        def ddegrad(s, c, t):
            Xslag = 0.0
            Plag = 0.0
            Xllag = 0.0

            if (t>c[7] and t<plyt_added): # if t > T for lysogenic phage
                Xllag = p.pastvalue(4, t - c[7], 0)

            if (t>(c[7]+plyt_added)): # if t > T for both phages
                Xslag = p.pastvalue(1,t-c[7],0)
                Plag = p.pastvalue(3,t-c[7],0)
                Xllag = p.pastvalue(4,t-c[7],0)

            if (s[1] < 1.0e-15):
                # if concentration of Xs drops below 1.0e-15 add time to the array
                Xs_extinction_times.append(t)

            g = array([0.0,0.0,0.0,0.0,0.0,0.0])

            # s[0] = S(t), s[1] = Xs(t), s[2] = Xi(t), s[3] = P(t), s[4] = Xl(t), s[5] = Pt(t)
            # S = - Xs*u*S/(Km+S)*(1/Y) - Xi*u*S/(Km+S)*(1/Y) - Xl*u*S/(Km+S)*(1/Y)
            g[0] = - s[1]*(c[0]*s[0])/(c[5]+s[0])*(1/c[6]) - s[2]*(c[0]*s[0])/(c[5]+s[0])*(1/c[6]) - s[4]*(c[0]*s[0])/(c[5]+s[0])*(1/c[6])
            # Xs = Xs*u*S/(Km+S)*(1/Y) - Ki*Xs*P - Ki*Xs*Pt
            g[1] = s[1]*c[0]*s[0]/(c[5]+s[0]) - c[3]*s[1]*s[3] - c[3]*s[1]*s[5]
            # Xi = Ki*Xs*P - Ki*Xs(t-T)P(t-T)
            g[2] = 0
            # P = b*Ki*Xs(t-T)P(t-T) - Ki*Xs*P
            g[3] = 0
            if (t>plyt_added): # after plyt_added hours lytic phage is added
                # Xi = Ki*Xs*P - Ki*Xs(t-T)P(t-T)
                g[2] = c[3]*s[1]*s[3] - c[3]*Xslag*Plag
                # P = b*Ki*Xs(t-T)P(t-T) - Ki*Xs*P - Ki*Xl*P - Ki*Xi*P
                g[3] = c[4]*c[3]*Xslag*Plag - c[3]*s[1]*s[3] - c[3]*s[4]*s[3] - c[3]*s[2]*s[3]
            # Xl = Xl*u*S/(Km+S) + Ki*Xs*Pt - q*Xl(t-T)
            g[4] = s[4]*(c[0]*s[0])/(c[5]+s[0]) + c[3]*s[1]*s[5] -c[8]*Xllag
            # Pt = b*q*Xl(t-T)  - Ki*Xs*Pt - Ki*Xl*Pt - Ki*Xi*Pt
            g[5] = c[4]*c[8]*Xllag - c[3]*s[1]*s[5] - c[3]*s[4]*s[5] - c[3]*s[2]*s[5]
            return g
        dde_camp = p.dde()
        dde_camp2 = p.dde()

        Xs_extinction_times = [] # stores times once Xs goes below 1.0e-15

        # Definte function to store history variables: state and gradient
        def ddesthist(g, s, c, t):
            return (s, g)
        global Xs0, Xi, P0, Xl, Pt0
        # Setting constants
        ddecons = array([u,S0,D,Ki,b,Km,Y,T,q])
        # Setting initial conditions S, Xs, Xi, P, Xl, Pt
        S = ddecons[1]
        ddeist = array([S0, Xs0, Xi, P0, Xl, Pt0])
        # Setting a state-scaling array for use in error control when values are very close to 0
        ddestsc = array([0,0,0,0,0,0])

        # Time: 0-plyt_added (Lysogenic phage only)
        dde_camp.initproblem(no_vars=6, no_cons=9, nlag=1, nsw=0, t0=0.0, t1=plyt_added+0.1, initstate=ddeist, c=ddecons, otimes= arange(0.0, plyt_added+0.1, 0.1), grad=ddegrad, storehistory=ddesthist)
        dde_camp.initsolver(tol=0.000005, hbsize=10000, dt=1.0, statescale=ddestsc)
        dde_camp.solve()

        # Time: from plyt_added to sim_length (Lysogenic + Lytic phage)
        # Setting initial conditions for the second simulation

        S_endval = dde_camp.data[:, 1][-1]
        Xs0_endval = dde_camp.data[:, 2][-1]
        Xi_endval = dde_camp.data[:, 3][-1]
        P0_new = 1.0e4 # introduction of P lytic
        Xl_endval = dde_camp.data[:, 5][-1]
        Pt0_endval = dde_camp.data[:, 6][-1]
        # Setting initial conditions S, Xs, Xi, P, Xl, Pt
        ddeist = array([S_endval, Xs0_endval, Xi_endval, P0_new, Xl_endval, Pt0_endval])
        # Setting a state-scaling array for use in error control when values are very close to 0
        ddestsc = array([0,0,0,0,0,0])

        dde_camp2.initproblem(no_vars=6, no_cons=9, nlag=1, nsw=0, t0=plyt_added, t1=sim_length, initstate=ddeist, c=ddecons, otimes= arange(plyt_added, sim_length, 0.1), grad=ddegrad, storehistory=ddesthist)
        dde_camp2.initsolver(tol=0.000005, hbsize=10000, dt=1.0, statescale=ddestsc)
        dde_camp2.solve()

        # Print values at the joining point between two simulations
        #print('At t= ' + str(dde_camp.data[:, 0][-1]) + ', Xs =' + str(dde_camp.data[:, 2][-1]))
        #print('At t= ' + str(dde_camp2.data[:, 0][0]) + ', Xs =' +str(dde_camp2.data[:, 2][0]))
        Xs_extinct = Xs_extinction_times[0]
        print('Xs went extinct at t= ' + str(Xs_extinct))

        # Calculation of average concentration
        def get_avg(values):
            less = values > 1.0e-15 # select only values above threshold
            return values[less].sum()/len(values[less])

        # Calculation of r value (ratio of xi/xs)
        xs_avg = get_avg(concatenate((dde_camp.data[:, 2],dde_camp2.data[:, 2])))
        #print('Average Xs = ' + str(xs_avg))
        xi_avg = get_avg(concatenate((dde_camp.data[:, 3],dde_camp2.data[:, 3])))
        #print('Average Xi = ' + str(xi_avg))
        xixs_ratio = round(xi_avg/xs_avg, 4)
        print('Xi/Xs ratio is ' + str(xixs_ratio))

        # Plot Test figure
        plt.style.use('ggplot') # set the global style
        xs, = plt.plot(concatenate((dde_camp.data[:, 0], dde_camp2.data[:, 0])), concatenate((dde_camp.data[:, 2],dde_camp2.data[:, 2])),  label=r'$X_S$')
        xi, = plt.plot(concatenate((dde_camp.data[:, 0], dde_camp2.data[:, 0])), concatenate((dde_camp.data[:, 3],dde_camp2.data[:, 3])),  label=r'$X_I$')
        plyt, = plt.plot(dde_camp2.data[:, 0], dde_camp2.data[:, 4], "r:",  label=r'$P$') # no need to concatenate since previous values = 0
        xl, = plt.plot(concatenate((dde_camp.data[:, 0], dde_camp2.data[:, 0])), concatenate((dde_camp.data[:, 5],dde_camp2.data[:, 5])),  label=r'$X_L$')
        pt, = plt.plot(concatenate((dde_camp.data[:, 0], dde_camp2.data[:, 0])), concatenate((dde_camp.data[:, 6],dde_camp2.data[:, 6])), "r--", label=r'$P_T$')
        f_size = 15 # set font size for plot labels
        plt.xlabel('Time (hours)', fontsize=f_size)
        plt.ylabel('Log concentration (particles/ml)', fontsize=f_size)
        plt.yscale('log')
        plt.axis([0,sim_length,1.0e-4,1.0e10])
        #plt.text(sim_length*0.2,8.0e8,'$P(t)$= '+str(plyt_added)+' h', fontsize=f_size) # display parameters
        plt.text(Xs_extinct,1.5e10,'$t=$ ' + str(round(Xs_extinct,3)), fontsize=f_size-1) # display parameters
        plt.text(0-0.5,1.5e10,'$r=$ ' + str(xixs_ratio), fontsize=f_size-1) # display parameters
        plt.tick_params(axis='both', labelsize=f_size)
        plt.vlines(Xs_extinct, 0,1.0e10, linewidth=0.5)
        # Plot substrate on the second y axis on top of the preivous figure
        plt2 = plt.twinx()
        plt2.grid(False)
        s, = plt.plot(concatenate((dde_camp.data[:, 0], dde_camp2.data[:, 0])), concatenate((dde_camp.data[:, 1],dde_camp2.data[:, 1])), 'black',  label=r'$S$')
        plt2.set_ylabel(r'Substrate (${\mu}$g/ml)', fontsize=f_size)
        plt2.set_yticks(linspace(0,S0, 3))
        plt2.tick_params(axis='both', labelsize=f_size)
        # Join legends from two separate plots into one
        plots = [xs,xi,plyt,xl,pt,s]
        plt.legend(plots, [plots_.get_label() for plots_ in plots],loc='best', fontsize= 'small', prop={'size': f_size})
        plt.tight_layout()
        plt.show()

        result[percentage]=[Xs_extinct,xixs_ratio]
    return result

parameter = 'u'
data = dde_sa(parameter, 95.0, 105.0, 10)
percentages = list(data.keys())
extinction_times = array(list(data.values()))[:,0]
r_values = array(list(data.values()))[:,1]
print(percentages)

# Plot figures

plt.style.use('ggplot') # set the global style
Xs_ext_plot, = plt.plot(percentages,extinction_times, label=r'$time$')
#xi, = plt.plot(concatenate((dde_camp.data[:, 0], dde_camp2.data[:, 0])), concatenate((dde_camp.data[:, 3],dde_camp2.data[:, 3])),  label=r'$X_I$')

f_size = 15 # set font size for plot labels
plt.xlabel('Percentage change in $'+parameter+'$', fontsize=f_size)
plt.ylabel('Time of $X_S$ extinction', fontsize=f_size)
#plt.yscale('log')
#plt.axis([0,20,1.0e-4,1.0e10])
#plt.text(sim_length*0.2,8.0e8,'$P(t)$= '+str(plyt_added)+' h', fontsize=f_size) # display parameters
#plt.text(Xs_extinct,1.5e10,'$t=$ ' + str(round(Xs_extinct,3)), fontsize=f_size-1) # display parameters
#plt.text(0-0.5,1.5e10,'$r=$ ' + str(xixs_ratio), fontsize=f_size-1) # display parameters
plt.tick_params(axis='both', labelsize=f_size)
fmt = '%.0f%%' # Format you want the ticks, e.g. '40%'
xticks = mtick.FormatStrFormatter(fmt)
ax = plt.gca()
ax.xaxis.set_major_formatter(xticks)


# Plot substrate on the second y axis on top of the preivous figure
plt2 = plt.twinx()
plt2.grid(False)
r_values_plot, = plt.plot(percentages, r_values, 'black',  label=r'$r$-value')
plt2.set_ylabel(r'$r$-value', fontsize=f_size)
plt2.set_yticks(linspace(min(r_values),max(r_values), 3))
plt2.tick_params(axis='both', labelsize=f_size)


# Join legends from two separate plots into one
p = [Xs_ext_plot,r_values_plot]
plt.legend(p, [p_.get_label() for p_ in p],loc='best', fontsize= 'small', prop={'size': f_size})
plt.tight_layout()
plt.show()
