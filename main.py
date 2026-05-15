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
mu = 3.178*10^-5   # Dynamic viscosity
V_gust = 80        # wind gust (150 kt)
T_atm = 226 

# Fin variables
c_t = 0.35          # tip chord   
c_r = 0.6        # root chord
b = 0.36           # fin span
t = 0.015          # thickness
LE_angle = 10     # LE angle
sweep = 37.5       # LE sweep
m_f = 1.5          # fin mass
cant_angle = 2     # difference of fin angle from vertical
dCL_da = 1.4792  

