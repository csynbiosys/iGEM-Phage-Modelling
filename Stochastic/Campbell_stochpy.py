import stochpy
from scipy import *
from numpy import *
import matplotlib.pyplot as plt

smod = stochpy.SSA(model_file='campbell_model.psc', dir='data')
smod.DoStochSim(end=10.0, mode='time')
smod.PlotSpeciesTimeSeries()
stochpy.plt.yscale('log')