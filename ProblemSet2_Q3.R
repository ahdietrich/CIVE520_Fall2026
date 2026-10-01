##Question 3 constants
rho_w <- 1000        # density of water, kg/m3
D_bar <- 0.9         # average drop diameter, mm
N_0 <- 8000          # 1/(m3 mm)
a <- 1500            # terminal velocity factor, 1/s
w <- 10              # updraft velocity, m/s

##Marshall-Palmer drop size distribution
#N(D) = N_0 * exp(-c*D), where c = 1/D_bar
c_MP <- 1/D_bar      # 1/mm
N_D <- function(D){N_0 * exp(-c_MP * D)}

#terminal velocity, linear form: v_t = a*D (a = alpha)  (D in mm -> m, gives m/s)
v_t <- function(D){a * D/1000}

#unit conversion: mm3 -> m3
mm3_m3 <- 1e-9

##(a) total number of drops
#N_T = integral from 0 to Inf of N(D) dD
N_T <- integrate(N_D, 0, Inf)$value          # drops per m3

##(b) liquid water content
#LWC = integral from 0 to Inf of (pi*D^3)/6 * rho_w * N(D) dD
#(pi*D^3)/6 is drop volume in mm3 -> convert to m3 so rho_w gives kg
LWC_kg <- integrate(function(D){(pi * D^3)/6 * mm3_m3 * rho_w * N_D(D)}, 0, Inf)$value  # kg/m3
LWC <- LWC_kg * 1000                          # g/m3

##(c) minimum drop size that falls out of the cloud
#drop falls when v_t > w  ->  a*D > w  ->  D_min = w/a
D_min <- (w/a) * 1000                         # mm
frac_D_min <- exp(-c_MP * D_min)              # fraction of drops larger than D_min

##(d) precipitation rate at the surface
#R = integral from D_min to Inf of (pi*D^3)/6 * fall speed * N(D) dD
#units: mm3 * m/s * 1/(m3 mm) * mm = mm3/(m2 s)
mmhr_conv <- mm3_m3 * 1000 * 3600             # mm3/(m2 s) -> m/s -> mm/s -> mm/hr

#net fall speed out of the cloud base (v_t - w)
R_net <- integrate(function(D){(pi * D^3)/6 * (v_t(D) - w) * N_D(D)}, D_min, Inf)$value * mmhr_conv  # mm/hr
R_net_in <- R_net/25.4                        # in/hr
