import numpy as np

# Gas constants (air)
gamma = 1.4 # specific heat ratio 
R = 287     # specific gas constant (J/kgK)

# Rocket variables
m_n = 3.69            # nosecone mass
l_n = 1.25         # length of nosecone
m_r = 275 - m_n     # rocket mass (without fins or nosecone)
l_r = 7.87 - l_n    # rocket length (without fins or nosecone
r = 0.125          # rocket radius
V = 924            # Rocket speed at qmax
rho = 0.3652       # Density at qmax
M1 = 2.8           # M at qmax     
p1 = 0.3398*100000 # Pressure at qmax
alpha = 10         # max AoA
mu = 3.178 * 10**-5   # Dynamic viscosity
V_gust = 80        # wind gust (150 kt)
T_atm = 226 

# Fin variables
c_t = 0.35          # tip chord   
c_r = 0.6        # root chord
b = 0.36           # fin span
t = 0.015          # thickness
LE_angle = 10     # LE angle (from wind direction, fin looks like a wedge)
sweep = 37.5       # LE sweep
m_f = 1.5          # fin mass
cant_angle = 2     # difference of fin angle from vertical
dCL_da = 1.4792  

qmax = 0.5 * rho * V**2 # max dynamic pressure
a = np.sqrt(gamma * R * T_atm) # speed of sound

# rocket subject to non-negligible wind disturbance
alpha_induced = np.atan(V_gust/V) # induced AoA


# for force normalisation, calculate planform area of fin
area_fin = 0.5 * (c_r + c_t) * b
mac = (2/3) * ((c_r**2 + c_r*c_t + c_t**2)/(c_r + c_t)) # mean aerodynamic chord (if trapezium became triangle)
M1_angle = np.rad2deg(np.asin(1/M1))

Re = (rho * V * mac)/mu

# Ackeret's theory (can only be applied for small angles but good first approx)
C_l = 4*np.deg2rad(alpha)/(np.sqrt(M1**2 - 1)) # lift coefficient approximation
print(C_l)

