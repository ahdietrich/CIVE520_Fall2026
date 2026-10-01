##defined constants
L_v = 2.5e6 #J/kg
L_s = 2.85e6 #J/kg
C_pw = 4216 #J/kgK
C_pd = 1004 #J/kgK
R_d = 287 #J/kgK
R_v = 461 #J/kgK
e_s0 = 611 #Pa
LvRv = L_v/R_v #K 
RvLv = R_v/L_v #1/K
epsilon = 0.622
g = 9.81 #m/s2

##question constants
T_0 <- 25 + 273.15   # surface temperature, K
q_0 <- 0.010         # parcel specific humidity, kg/kg
p_0 <- 101300        # surface pressure, Pa
G_a <- 7e-3          # ambient lapse rate, K/m
G_d <- g / C_pd        # dry adiabatic lapse rate, K/m


##define functions
#Saturation vapor pressure
e_s <- function(T_parcel){e_s0 * exp(LvRv * (1/273.15 - 1/T_parcel))}

#saturation specific humidiity 
q_s <- function(T_parcel,P_parcel) {epsilon * e_s(T_parcel)/P_parcel}

#Temperature and pressure profiles
T_env <- function(z){T_0 - G_a*z}
pz <- function(z){p_0 * (T_env(z)/T_0)^(g/(R_d * G_a))}

#Moist adibatic lapse rate (from hint)
G_m <- function(T_parcel,P_parcel){
  dqsdT <- L_v * q_s(T_parcel,P_parcel)/(R_v * T_parcel^2)
  g/(C_pd + L_v * dqsdT)
}

#trapezoid rule to integrate CAPE
trapz <- function(x,y){sum(diff(x) * (head(y,-1) + tail(y,-1))/2)}


#Lift parcel from 1:12000m
dz <- 1
z <- seq(0,12000,by=dz)
n <- length(z)
T_p <- numeric(n); T_p[1] <- T_0
z_LCL <- NA
T_LCL <- NA

for (i in 2:n) {
  if (is.na(z_LCL)) {
    # unsaturated: dry adiabatic
    T_p[i] <- T_p[i - 1] - G_d * dz
    if (q_s(T_p[i], pz(z[i])) <= q_0) {
      z_LCL <- z[i]; T_LCL <- T_p[i]
    }
  } else {
    # saturated: moist adiabatic
    T_p[i] <- T_p[i - 1] - G_m(T_p[i - 1], pz(z[i - 1])) * dz
  }
}

T_e <- T_env(z)
buoy <- T_p - T_e

#saturated lapse rate at LCL
G_m_LCL <- G_m(T_LCL, pz(z_LCL))

#level of free convection (LFC)
z_LFC <- z[which(z > z_LCL & buoy >0)[1]]

##CAPE from LFC to 12000m
#trapezoid rule to estimate integral
m <- z >= z_LFC
CAPE <- trapz(z[m],g * buoy[m]/T_e[m])
updraft_thr <- sqrt(2 * CAPE)
updraft_act <- updraft_thr/2

##saturation q_S of parcel
q_s_p <- q_s(T_p, pz(z))


##plots (claude generated code)
##Q2c: environment vs parcel temperature
plot(T_e - 273.15, z/1000, type = "l", lwd = 2, col = "steelblue",
     xlim = range(c(T_e, T_p) - 273.15),
     xlab = "Temperature (\u00B0C)", ylab = "Height z (km)",
     main = "Q2c: Environment vs. parcel temperature")
lines(T_p - 273.15, z/1000, lwd = 2, col = "darkorange")
abline(h = z_LCL/1000, lty = 2, col = "gray40")
abline(h = z_LFC/1000, lty = 3, col = "gray40")
text(-55, z_LCL/1000 + 0.2, sprintf("LCL = %.2f km", z_LCL/1000), adj = 0)
text(-55, z_LFC/1000 + 0.2, sprintf("LFC = %.2f km", z_LFC/1000), adj = 0)
legend("topright", c("Environment", "Parcel"),
       col = c("steelblue", "darkorange"), lwd = 2)
grid()



plot(q_s_p * 1000, z/1000, type = "l", lwd = 2, col = "steelblue",
     xlab = expression("Saturation specific humidity " * q[s] * " (g/kg)"),
     ylab = "Height z (km)",
     main = "Q2d: Parcel saturation specific humidity")
abline(v = q_0 * 1000, lty = 2, col = "gray40")
abline(h = z_LCL/1000, lty = 3, col = "gray40")
legend("topright", c(expression(q[s] * " of parcel"), "q = 10 g/kg"),
       col = c("steelblue", "gray40"), lty = c(1, 2), lwd = c(2, 1))
grid()
