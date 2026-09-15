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
RvLv = R_v/L_v #1/K
epsilon = 0.622

##question constants
T_a = 303.15 # air temperature (K)
T_w = 298.15 # wet-bulb temperature (K)
P_sfc = 101325 #surface pressure (Pa)

#Clausius-Clapeyron Eqn
e_sT = e_s0*exp(LvRv*(1/T_0 - 1/T_a)) #4358.229 
e_sTw = e_s0*exp(LvRv*(1/T_0 - 1/T_w)) #3228.665

##Saturated specific humidty
q_s = epsilon*(e_sTw/P_sfc) #0.01981969
                
##Wet-bulb Temperature eqn rearrange to get specific humidty 
q = q_s - (C_pd*(T_a - T_w))/L_v #0.01781169

##vapor pressure
e = (q*P_sfc)/epsilon #2901.558

##Dewpoint Temperature
T_d = ((1/T_0)-RvLv*log(e/e_s0))^-1 #296.4092 

##Relative Humidity
RH = e/e_sT # 0.6657655

##air density
rho_v = e/(R_v*T_a) #0.02076217
rho_d = (P_sfc - e)/(R_d*T_a) #1.131251

rho = rho_v + rho_d  #1.152013



