##defined constants
L_v = 2.5e6 #J/kg
L_s = 2.85e6 #J/kg
C_pw = 4216 #J/kgK
C_pd = 1004 #J/kgK
R_d = 287 #J/kgK
R_v = 461 #J/kgK
e_s0 = 611 #Pa
T_0 = 273.15 #K
LvRv = L_v/R_v #K 
LsRv = L_s/R_v #K
epsilon = 0.622

##question constants
T_sfc = 290 #K
P_sfc = 101325 #Pa
T_10km = 225 #K
P_10km = 26500 #Pa

#Clausius-Clapeyron Eqn
e_sfc = e_s0*exp(LvRv*(1/T_0 - 1/T_sfc)) #1936.536
e_10km = e_s0*exp(LsRv*(1/T_0 - 1/T_10km)) #4.814481


#vapor density from ideal gas law
rho_sfc = e_sfc/(R_v*T_sfc) #0.01448527
rho_10km = e_10km/(R_v*T_10km) #4.641582e-05



