# Reactions
######################### S
S_consumption_by_Xs:
    S > $pool
    Xs*u*S/(Km+S)*(1/Y)

S_consumption_by_Xi:
    S > $pool
    Xi*u*S/(Km+S)*(1/Y)

S_consumption_by_Xl:
    S > $pool
    Xl*u*S/(Km+S)*(1/Y)

######################### Xs
growth_Xs:
    Xs > Xs + Xs
    Xs*u*S/(Km+S)

death_Xs:
    Xs > $pool
    Xs*D

infection_by_P:
    Xs + P > Xi
    Ki*Xs*P

infection_by_Pt:
    Xs + Pt > Xl
    Ki*Xs*Pt

######################### P

P_absorbed_onto_Xl:
    P + Xl > Xl
    Ki*Xl*P

######################### Xl
growth_Xl:
    Xl > Xl + Xl
    Xl*u*S/(Km+S)

death_Xl:
    Xl > $pool
    Xl*D

######################### Pt
    
Pt_absorbed_onto_Xl:
    Xl + Pt > Xl
    Ki*Xs*Pt

######################### Xi
death_Xi:
    Xi > $pool
    Xi*D
    
lysis_of_Xi:
    Xi > {98}P
    Xi*(1/0.5) #Lysis rate is cells/hr


# Variable species
S = 10.0
Xi = 1000
Xs = 1.0e4
Xl = 0
P = 1.0e2
Pt = 1.0e4

# Parameters
u = 0.738 #(Levin et al., 1977)
#S0 = 100.0 #ug/ml(Levin et al., 1977)
D = 0.20 # h-1  HERE USED AS COMMON DEATH RATE FOR BACTERIA
Ki = 6.24e-8 #ml/h (Levin et al., 1977)
b = 98.0 #(Levin et al., 1977)
Km = 4.0 #4 ug/ml(Levin et al., 1977)
Y = 7.40e4 #(Levin et al., 1977)
#T = 0.5 #h-1 (Levin et al., 1977)
#Xs0 = 1.0e6 # cells/ml starting levels of cells (Levin et al., 1977)
#P0 = 0 # particles/ml starting levels of cells (Levin et al., 1977)
#q = 0.001 # induction rate (...) MIGHT BE REMOVED LATER AT ALL
