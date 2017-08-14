# Reactions
######################### S
S_inflow:
    $pool > S
    D*S0
    
S_outflow:
    S > $pool
    D*S

######################### Xs
growth_Xs:
    Xs + S > {2.0}Xs
    Xs*u*S/(Km+S)

death_Xs:
    Xs > $pool
    Xs*D

infection_by_P:
    Xs + P > Xi
    Ki*Xs*P


######################### P

outflow_P:
    P > $pool
    P*D


######################### Xi
death_Xi:
    Xi > $pool
    Xi*D

lysis_of_Xi:
    Xi > {98}P
    Xi*1/T #Lysis rate is cells/hr
    
S_cons_by_Xi:
    Xi + S > Xi
    Xi*u*S/(Km+S)


# Variable species
S = 3000.0 # ug/ml
Xi = 0
Xs = 1.0e4
P = 1.0e2


# Parameters
u = 0.738 #(Levin et al., 1977)
S0 = 100.0 #ug/ml(Levin et al., 1977)
D = 0.20 # h-1  HERE USED AS COMMON DEATH RATE FOR BACTERIA
Ki = 6.24e-8 #ml/h (Levin et al., 1977)
#b = 98.0 #(Levin et al., 1977)
Km = 4.0 #4 ug/ml(Levin et al., 1977)
Y = 7.40e4 #(Levin et al., 1977)
T = 0.5 #h-1 (Levin et al., 1977)
