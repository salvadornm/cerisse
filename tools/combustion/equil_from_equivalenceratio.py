"""
Adiabatic flame temperature including solid carbon formation
============================================================

Adiabatic flame temperature and equilibrium composition for a fuel/air mixture
as a function of equivalence ratio, including formation of solid carbon.

Requires: cantera >= 2.5.0, matplotlib >= 2.0

.. tags:: Python, equilibrium, combustion, multiphase
"""

import cantera as ct
import numpy as np
import sys
import csv

##############################################################################
# Edit these parameters to change the initial temperature, the pressure, and
# the phases in the mixture.

T = 298.0
P = 101325.0
phi = 0.5
print(" --------------------- ")
print(" Initial Conditions  Mixture")
print("    P= ",P, "[Pa] amd T = ",T, "[K]")
print("    Equivalence Ratio= ",phi)

# phases
gas = ct.Solution('Li.yaml')

# the phases that will be included in the calculation, and their initial moles
mix_phases = [(gas, 1.0)]

# gaseous fuel species
fuel_species = 'H2'

##############################################################################

mix = ct.Mixture(mix_phases)

# set the gas state
gas.set_equivalence_ratio(phi, fuel_species, 'O2:1.0, N2:3.76')

# create teh mixture
mix = ct.Mixture(mix_phases)
mix.T = T
mix.P = P

print(" Unburn Mixture ...")
print(mix.report())


# equilibrate the mixture adiabatically at constant P
mix.equilibrate('HP', solver='gibbs', max_steps=1000)

# some array to store species
xeq = np.zeros(mix.n_species)
xeq[:] = mix.species_moles    # equilbrium composition
Teq = mix.T                   # equilibrium T

print(" Burn Mixture ...")
print(mix.report())

    
