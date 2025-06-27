import cantera as ct
import numpy as np
import matplotlib.pyplot as plt
import csv
import sys
import os

# Path to directory containing nozzle_functions.py
sys.path.append("/Users/snm/codes/cerisse/tools/nozzle")

from nozzle_functions import nozzle_functions


# Free-Stream ---------

# 1. Define gas mixture and mechanism
gas = ct.Solution("mechanism.yaml")  # Use a mechanism that includes CO2/Ar/N2 combustion

# 2.  Free strram
T_oo= 227.0  # Initial temperature [K]
p_oo = 284.0  # Initial pressure [Pa]
Mach = 2.0

# 3. Define the mixture composition as mass fractions
mix = {"CO2": 0.96, "AR": 0.04}

# 4. Set the state freestream conditions
gas.TPY = T_oo, p_oo, mix

print(" Mixture ...")
print(gas.report())

# Compute speed of sound
c_oo = gas.sound_speed  # [m/s]
rho_oo = gas.density    
u_oo    = c_oo*Mach;  
kin_oo  = 0.5*rho_oo*u_oo*u_oo
eint_oo = gas.int_energy_mass  # [J/kg] internal energy

print(f"Speed of sound: {c_oo:.3f} m/s")

gamma = gas.cp_mass / gas.cv_mass
print(f"Adiabatic index (γ): {gamma:.4f}")

# Stagnation properties
P0 = nozzle_functions.Pstag(p_oo, Mach, gamma)
T0 = nozzle_functions.Tstag(T_oo, Mach, gamma)

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
print("Stagnation conditions:")
print("P0 =", P0)
print("T0 =", T0)
print(" ------------------------------------------------")

# SRP 
print("Stagnation conditions SRP:")
P0srp = P0
T0srp = 5.0*T0

# 1. Define gas mixture and mechanism
gas_srp = ct.Solution("mechanism.yaml")  # Use a mechanism that includes CO2/Ar/N2 combustion
mix_srp = {"N2": 1.0}
gas_srp.TPY = T0srp, P0srp, mix_srp

print("P0srp =", P0srp)
print("T0srp =", T0srp)
print(" ------------------------------------------------")
gamma_srp = gas_srp.cp_mass / gas_srp.cv_mass
print(f"Adiabatic index (γ): {gamma_srp:.4f}")
Pt =nozzle_functions.Pchok(P0srp, gamma_srp)
Tt = nozzle_functions.Tchok(T0srp, gamma_srp)
c_srp    = gas_srp.sound_speed
u_srp    = c_srp
print("Throat conditions: (SRP)")
print("Pt =", Pt)
print("Tt =", Tt)
print("c_srp =", c_srp)
print("u_srp =", u_srp)
print(" ------------------------------------------------")
# massflow = nozzle_functions.masschok(T0srp, P0srp, gamma_srp)
# print("Mass flow rate at throat (SRP):", massflow)
# Mach_srp = 2.  # Mach number at nozzle exit
# Pe = nozzle_functions.Pnozz(Mach_srp, P0srp, gamma_srp)
# print("SRP Nozzle exit pressure =", Pe)


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
#print(" // SRP nozzle exit pressure (based on Mach number", Mach_srp, ")")
#print(" static constexpr Real Pe_srp  = ", Pe, ";")
print(" //-----------------------------------------------------------------")










