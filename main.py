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
alpha_ind = np.rad2deg(np.atan(V_gust/V)) # induced AoA
alpha_eff = alpha + alpha_ind
# print(alpha_eff)

# for force normalisation, calculate planform area of fin
S_area = 0.5 * (c_r + c_t) * b # area of a trapezium
mac = (2/3) * ((c_r**2 + c_r*c_t + c_t**2)/(c_r + c_t)) # mean aerodynamic chord (if trapezium became triangle)
M1_angle = np.rad2deg(np.asin(1/M1))
A_r = b**2/S_area

Re = (rho * V * mac)/mu # Reynolds number at MAC

# Ackeret's theory (can only be applied for small angles but good first approx)
# C_l = 4*np.deg2rad(alpha)/(np.sqrt(M1**2 - 1)) # lift coefficient approximation
# print(C_l)

# M normal to LE
# print(M1 * np.cos(np.deg2rad(sweep))) 
# shock attaches cleanly to LE so 2D strip theory can be used

# Sectional C_L using shock expansion theory
# treat fins as flat plate turning flow by alpha_eff

# helper functions
# beta in radians
def beta2theta(beta, M, gamma): 
    num = M**2 * (np.sin(beta))**2 - 1
    denom = M**2 * (gamma + np.cos(2 * beta)) + 2
    theta = np.atan((2 / np.tan(beta)) * (num / denom))
    return theta

def OS_pressure_ratio(gamma, M, beta):
    num = 2 * gamma * (M**2 * (np.sin(beta))**2 - 1)
    denom = gamma + 1
    p_ratio = 1 + num / denom
    return p_ratio

def M_up2M_down(gamma, M, beta, theta):
    M1_norm = M * np.sin(beta)
    coeff = (gamma - 1) / 2
    num = 1 + coeff * M1_norm**2
    denom = gamma * M1_norm**2 - coeff
    M2_norm = np.sqrt(num / denom)
    M_down = M2_norm / np.sin(beta - theta)
    return M_down

def prandtl_meyer(gamma, M):
    term1 = np.sqrt((gamma + 1) / (gamma - 1))
    term2 = np.sqrt(M**2 - 1)
    term3 = np.atan(term1 * term2)
    nu_M = term1 * term3 - np.atan(term2)
    return nu_M

def isen_pressure_ratio(gamma, M1, M2):
    coeff = (gamma - 1) / 2
    num = 1 + coeff * M1**2
    denom = a + coeff * M2**2
    power = gamma / (gamma - 1)
    p_ratio = (num / denom)**power
    return p_ratio