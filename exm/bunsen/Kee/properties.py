import cantera as ct
import numpy as np
import matplotlib.pyplot as plt
import csv

# 1. Define gas mixture and mechanism
gas = ct.Solution("Kee.yaml")  # Use a mechanism that includes H2/Air combustion

# 2. Set  conditions
T = 300.0       # Initial temperature [K]
P = 101325.0   # Initial pressure [Pa]

# 3. Define the mixture composition as mass fractions
mix = {"H2": 0.875, "CH4": 0.125}

# 4. Set the state
gas.TPY = T, P, mix

print(" Mixture ...")
print(gas.report())

# Compute speed of sound
a = gas.sound_speed  # [m/s]

print(f"Speed of sound: {a:.3f} m/s")



