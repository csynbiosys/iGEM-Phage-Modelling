from numpy import *

# Importing data points (time, values) from Fig 3a from Krysiak-Baltyn et al. (2016)
bacteria = genfromtxt('data/3a_bateria_non-log.csv', delimiter=',')
phage = genfromtxt('data/3a_phage_non-log.csv', delimiter=',')

print(phage)