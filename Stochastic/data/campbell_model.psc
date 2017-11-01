# Reactions

######################### Xs
growth_Xs:
    Xs > Xs + Xs
    Xs*u*(1-Xs/C)

washout_Xs:
    Xs > $pool
    Xs*D

infection_by_P:
    Xs + P > Xi #Typo: Xs + P > Xs changed to Xs + P > Xi
    Xs*P*Ki #Rate of infection is dependent on Xs and P


lysis_of_Xi:
    Xi > {98}P
    Xi*(1/0.5) #Lysis rate is cells/hr

#washout_P:
#    P > $pool
#    P*D
    
#deactivation_of_P:
#    P > $pool
#    P*dp
    
######################### Xi
#washout_Xi:
#    Xi > $pool
 #   Xi*D



# Variable species
Xi = 0   # cells/ml
Xs = (1.0e2) * 5    # cells/ml
P = (1.0e3) * 5   # particles/ml

# Parameters
#Na = 6.02e23  # mol-1
V = 5 # ml
u = 0.738 # h-1
D = 0.20 # h-1
dp = 0.1 # h-1
Ki = 6.24e-8 # ml attacks per phage particle per cell per hour (Levin et al., 1977)
b = 98.0 # (Levin et al., 1977)
C = 7.0e9 # cells
#T = 0.5 #h-1 (Levin et al., 1977)
