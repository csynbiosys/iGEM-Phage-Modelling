import stochpy
from scipy import *
from numpy import *
import matplotlib.pyplot as plt

smod = stochpy.SSA(model_file='bact_phage_model.psc', dir='data')
smod.DoStochSim(end=500.0, mode='time')
smod.PlotSpeciesTimeSeries()
stochpy.plt.yscale('log')

#S = smod.data_stochsim
#t= reshape(S.time, len(S.time))
#print(S.species[:,0])

#plt.plot(t, S.species[:,1], label='1')
#plt.legend()
#plt.xlabel('Time (hours)')
#plt.ylabel('Species')
##plt.yscale('log')
#plt.axis([-5,250,1.0e-4,1.0e10])
#plt.tick_params( axis='both',labelsize=12)
#plt.show()
