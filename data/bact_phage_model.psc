# Reactions

S_consumption_by_Xs:
    S > $pool
    Xs*u*S/(Km+S)
    
S_consumption_by_Xi:
    S > $pool
    Xi*u*S/(Km+S)
    

growth_Xs:
    Xs > Xs + Xs
    Xs*u*S/(Km+S)
    
infection_Xi:
    Xs > Xi
    Ki*Xs*P
    
death_Xi:
    Xi > $pool
    Xi*D
    
lysis_P:
    P > P + P
    P*b*Ki*Xs
    
infection_P:
    P > $pool
    Ki*Xs*P
    

    

    
    



# Variable species
S = 30
Xi = 0 # no lytic bacteria present at the start
Xs = 1.0e4
P = 1.0e4
#Pt = 1.0e6 # particles/ml of temperate phage
#Xl = 0 # no lysogenic bacteria present at the start

# Parameters
u = 0.738 #(Levin et al., 1977)
S0 = 100.0 #ug/ml(Levin et al., 1977)
D = 0.20 # h-1  HERE USED AS COMMON DEATH RATE FOR BACTERIA
Ki = 6.24e-8 #ml/h (Levin et al., 1977)
b = 98.0 #(Levin et al., 1977)
Km = 4.0 #4 ug/ml(Levin et al., 1977)
Y = 7.40e4 #(Levin et al., 1977)
T = plyt_added + 0.5 #h-1 (Levin et al., 1977)
Xs0 = 1.0e6 # cells/ml starting levels of cells (Levin et al., 1977)
P0 = 0 # particles/ml starting levels of cells (Levin et al., 1977)
q = 0.001 # induction rate (...) MIGHT BE REMOVED LATER AT ALL
