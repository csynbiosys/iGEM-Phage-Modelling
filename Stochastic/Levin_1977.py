import stochpy
from scipy import *
from numpy import *
import matplotlib.pyplot as plt

smod = stochpy.SSA(model_file='levin_1977.psc', dir='data')
smod.DoStochSim(end=0.1, mode='time', IsTrackPropensities=True)
smod.PlotSpeciesTimeSeries() # plot time series of species
stochpy.plt.yscale('log')
smod.PlotPropensitiesTimeSeries() # plot time series of propensities
#smod.PlotWaitingtimesDistributions()


