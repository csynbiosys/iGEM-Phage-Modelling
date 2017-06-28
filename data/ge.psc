# Reactions

transcription:
    D > M + D
    u

degradation_m:
    M > $pool
    dm*M

translation:
    M > M + P
    v*M

degradation_p:
    P > $pool
    dp*P

# Fixed species

# Variable species
D= 1
M= 0
P= 0

# Parameters
u= 0.01
dm= 0.00385  # log(2)/180
v= 0.04  
dp= 0.00019 # log(2)/3600
