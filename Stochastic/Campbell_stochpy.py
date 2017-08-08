import stochpy
from scipy import *
from numpy import *
import matplotlib.pyplot as plt

smod = stochpy.SSA(model_file='campbell_model_bradley.psc', dir='data')
smod.DoStochSim(quiet=False, end=2.0, mode='time')
smod.PlotSpeciesTimeSeries()
stochpy.plt.yscale('log')