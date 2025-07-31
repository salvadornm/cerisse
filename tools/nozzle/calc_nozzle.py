import numpy as np
import math
import os
import sys

# Path to directory containing nozzle_functions.py
sys.path.append("/Users/snm/codes/cerisse/tools/nozzle")

from nozzle_functions import nozzle_functions

# gas
gamma = 1.4
Rgas  = 287.0  # [J/(kg*K)] specific gas constant for air
Cv   = Rgas/(gamma-1.0)  # [J/(kg*K)] specific heat at constant volume

# Free-stream conditions
p_oo    = 574.56;   # [Pa] free-stream pressure   
T_oo    = 65.0;     # [K]  free-stream temperature 
Mach    = 4.6       # Mach number

# derived free-stream conditions
rho_oo  = p_oo/(Rgas*T_oo)
c_oo    = np.sqrt(gamma*Rgas*T_oo)
u_oo    = c_oo*Mach;  
eint_oo = Cv*T_oo
kin_oo  = 0.5*rho_oo*u_oo*u_oo

print(" -----------  Free-stream conditions  -----------")
print("Mach =", Mach)
print("P_oo =", p_oo)
print("T_oo =", T_oo)
print("rho_oo =", rho_oo)
print("c_oo =", c_oo)
print("u_oo =", u_oo)  
print("eint_oo =", eint_oo)
print("kin_oo =", kin_oo) 
print(" ------------------------------------------------")
# stagnation pressure and temperature
P0 = nozzle_functions.Pstag(p_oo, Mach, gamma)
T0 = nozzle_functions.Tstag(T_oo, Mach, gamma)
print("Stagnation conditions:")
print("P0 =", P0)
print("T0 =", T0)
print(" ------------------------------------------------")

# SRP stagation P and T (relative to free-stream conditions)
P0srp = 7724.0*p_oo
T0srp = 5.34*T_oo

print("Stagnation conditions SRP:")
print("P0srp =", P0srp)
print("T0srp =", T0srp)
print(" ------------------------------------------------")

# compute conditions at throat
Pt =nozzle_functions.Pchok(P0srp, gamma)
Tt = nozzle_functions.Tchok(T0srp, gamma)
c_srp    = np.sqrt(gamma*Rgas*Tt)
u_srp    = c_srp
print("Throat conditions: (SRP)")
print("Pt =", Pt)
print("Tt =", Tt)
print("c_srp =", c_srp)
print("u_srp =", u_srp)
print(" ------------------------------------------------")
massflow = nozzle_functions.masschok(T0srp, P0srp, gamma)
print("Mass flow rate at throat (SRP):", massflow)

Mach_srp = 2.94  # Mach number at nozzle exit
Pe = nozzle_functions.Pnozz(Mach_srp, P0srp, gamma)
print("SRP Nozzle exit pressure =", Pe)


print(" for cerisse prob.h file (copy the  C++ lines below):")
print(" //------ values obtained from calc_nozzle.py  (SI units)---------")
print(" // based on Mach number", Mach, "and gamma", gamma)
print(" // and free-stream conditions P_oo =", p_oo, "[Pa] and T_oo =", T_oo, "[K]")
print(" // Staganation conditions")
print(" static constexpr Real P0  = ", P0, ";")
print(" static constexpr Real T0  = ", T0, ";")
print(" // Free-stream conditions")
print(" static constexpr Real p_oo  = ", p_oo, ";")
print(" static constexpr Real T_oo  = ", T_oo, ";")
print(" static constexpr Real rho_oo  = ", rho_oo, ";")
print(" static constexpr Real c_oo  = ", c_oo, ";")
print(" static constexpr Real u_oo  = ", u_oo, ";")
print(" static constexpr Real eint_oo  = ", eint_oo, ";")
print(" static constexpr Real kin_oo  = ", kin_oo, ";")
print(" // SRP Staganation conditions")
print(" static constexpr Real P0srp  = ", P0srp, ";")
print(" static constexpr Real T0srp  = ", T0srp, ";")
print(" // SRP throat conditions (choked flow)")
print(" static constexpr Real Pt  = ", Pt, ";")
print(" static constexpr Real Tt  = ", Tt, ";")
print(" static constexpr Real u_srp  = ", u_srp, ";")
print(" // SRP nozzle exit pressure (based on Mach number", Mach_srp, ")")
print(" static constexpr Real Pe_srp  = ", Pe, ";")
print(" //-----------------------------------------------------------------")
