"""
Computes sectional and 3D lift/drag coefficients via shock-expansion theory and
linear-slope methods, plus steady-state roll rate, centripetal load, and
leading-edge stagnation temperature. 
 
References:
    [1] Ribera Vicent, M. Compressible Aerodynamics, AERO50001 notes, Imperial.
    [2] Anderson, J.D. Fundamentals of Aerodynamics, 6th ed., McGraw-Hill.
    [3] Niskanen, S. OpenRocket Technical Documentation v13.05 (2013), §3.3.
    [4] Gudmundsson, S. General Aviation Aircraft Design, 2nd ed., 2022, Ch. 15.
    [5] Schlichting, H. Boundary-Layer Theory, 8th ed.
"""

import numpy as np
from scipy.optimize import brentq, minimize_scalar
from scipy.integrate import quad

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

# Derived flight conditions
qmax = 0.5 * rho * V**2 # max dynamic pressure
a = np.sqrt(gamma * R * T_atm) # speed of sound

# rocket subject to non-negligible wind disturbance
alpha_ind = np.rad2deg(np.atan(V_gust/V)) # induced AoA
alpha_eff = alpha + alpha_ind
# print(alpha_eff)

# for force normalisation, calculate planform area of fin
S_area = 0.5 * (c_r + c_t) * b # area of a trapezium
mac = (2/3) * ((c_r**2 + c_r*c_t + c_t**2)/(c_r + c_t)) # mean aerodynamic chord (if trapezium became rectangle)
mu_1 = np.asin(1/M1)  # Mach angle in radians
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
    term3 = np.atan(term2 / term1)
    nu_M = term1 * term3 - np.atan(term2)
    return nu_M # radians

def isen_pressure_ratio(gamma, M1, M2):
    coeff = (gamma - 1) / 2
    num = 1 + coeff * M1**2
    denom = 1 + coeff * M2**2
    power = gamma / (gamma - 1)
    p_ratio = (num / denom)**power
    return p_ratio

# since shock wave angle is difficult to isolate from the beta-theta relation,
# i will use scipy library for a numerical root finding algorithm 
# since effective AoA > LE_angle, lower front will have an oblique shock 
# (flow turning into itself), upper front PM expansion fan (flow turning away)

# lower surface
theta1 = np.deg2rad(LE_angle + alpha_eff)

# beta2theta peaks at some beta_max then falls back to 0 at pi/2,
# so both endpoints of [mu_1, pi/2] have the same sign — bracket won't work.
# Find the beta of maximum deflection first, then bracket the weak-shock root.
beta_max = minimize_scalar(lambda b: -beta2theta(b, M1, gamma),
                           bounds=(mu_1 + 1e-6, np.pi/2 - 1e-6),
                           method='bounded').x
beta = brentq(lambda b: beta2theta(b, M1, gamma) - theta1,
              mu_1 + 1e-6, beta_max - 1e-6)

# print(np.rad2deg(beta), np.rad2deg(beta_max))

pRatio_lower = OS_pressure_ratio(gamma, M1, beta)
M2_lower = M_up2M_down(gamma, M1, beta, theta1)

# upper surface
theta2 = np.deg2rad(alpha_eff - LE_angle)
nu1 = prandtl_meyer(gamma, M1)
nu_up = nu1 + theta2

# print(np.rad2deg(nu1), np.rad2deg(nu_up))

M2 = brentq(lambda M: prandtl_meyer(gamma, M) - nu_up, 1.001, 25.0)
pRatio_upper = isen_pressure_ratio(gamma, M1, M2)

# print(M2, pRatio_upper, pRatio_lower)

# from shock-expansion theory
Cp_lower = (2 / (gamma * M1**2)) * (pRatio_lower - 1)
Cp_upper = (2 / (gamma * M1**2)) * (pRatio_upper - 1)

Cl_sec = (Cp_lower - Cp_upper) * np.cos(np.deg2rad(alpha_eff)) # assuming infinite span, conservative upper bound
Cd_lift_sec = (Cp_lower - Cp_upper) * np.sin(np.deg2rad(alpha_eff)) # lift induced wave drag 2D
Cl_3D = dCL_da * np.deg2rad(alpha_eff) # finite-span C_L, closer to value at tip

print(f"C_L (2D shock-exp, upper-bound) = {Cl_sec:.4f}")
print(f"C_L (3D linear slope, lower-bound) = {Cl_3D:.4f}")

# wave drag + skin friction
tau = t / mac
Cd_thick = 4 * tau**2 / np.sqrt(M1**2 - 1) # thin aerofoil wave drag for fin thickness

Cf_emp = 0.0592 / Re**0.2 # Prandtl-Schlichting 1/7 power law for flat plate for turbulent, incompressible BL
Cf_comp = Cf_emp * (1 + 0.144 * M1**2)**(-0.65) # correction factor for compressiblity from Frankl-Voishel empirical correction

Cd_friction = 2 * Cf_comp # 2 because both sides of fin are wetted

Cd_total = Cd_lift_sec + Cd_thick + Cd_friction

print(f"C_D total = {Cd_total:.4f} "
      f"(lift-induced wave = {Cd_lift_sec:.4f}, "
      f"thickness wave = {Cd_thick:.4f}, friction = {Cd_friction:.4f})")

# Centripetal force calculation
chord = lambda y: c_r + (c_t - c_r) * y / b   # local chord
I1 = quad(lambda y: chord(y) * (r + y), 0, b)[0]  # first moment (forcing arm)
I2 = quad(lambda y: chord(y) * (r + y)**2, 0, b)[0]  # second moment (damping arm)
y_eff = I2 / I1
p_ss = V * np.tan(np.deg2rad(cant_angle)) / y_eff

# spanwise mass centroid of trapezium
y_bar = (b / 3) * (c_r + 2 * c_t) / (c_r + c_t)
r_cg  = r + y_bar

F_centripetal = m_f * p_ss**2 * r_cg # F = m * w^2 * r

print(f"Steady-state roll rate = {p_ss:.2f} rad/s ({p_ss/(2*np.pi):.2f} rev/s)")
print(f"Centripetal force (mass centroid) = {F_centripetal:.1f} N")

# Temperature at LE
T_stagnation = T_atm * (1 + 0.5*(gamma-1) * M1**2) # isentropic relation

print(f"Stagnation temperature at LE = {T_stagnation:.1f} K "
      f"({T_stagnation-273.15:.0f} deg C)")