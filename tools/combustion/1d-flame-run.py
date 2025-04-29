###############################################################
# ADIABATIC_FLAME - A freely-propagating, premixed flat flame 
#                   of fuel in AIR (O2 + N2) 
###############################################################

#import :
import cantera as ct
from cantera import *
import numpy as np
import csv
import matplotlib.pyplot as plt

#################################################################
# Prepare your run
#################################################################
# Parameter values :

# Output format
specformat = 'mass'   # mole or mass
units      = 'MKS'    #MKS or CGS

# 1. Define gas mixture and mechanism
gas = ct.Solution("Li.yaml")  # Use a mechanism that includes H2/Air combustion

# Define fuel and oxidizer mole fractions for a stoichiometric mixture
fuel = "H2"
oxidizer = {"O2": 1.0, "N2": 3.76}  # Air composition (O2/N2 ratio of 1:3.76)

# 2. Set initial conditions
P0           = ct.one_atm       # 1 atm [Pa]
p            = P0               # pressure [Pa]
tin          = 298.0            # unburned gas temperature [K]
phi          = 0.5              # Eq. ratio [-] 

# Refined grid at inlet and outlet, 6 points in x-direction :
domain_size  = 0.04             # Domain size [m] 
#initial_grid = domain_size*np.array([0.0, 0.001, 0.01, 0.02, 0.029, 0.03],'d')/0.03 # m

initial_grid = domain_size*np.array([0.0, 0.001, 0.3, 0.66, 0.98, 1.0]) # m


#Set tolerance properties
tol_ss    = [1.0e-8, 1.0e-9]        # [rtol atol] for steady-state problem
tol_ts    = [1.0e-8, 1.0e-9]        # [rtol atol] for time stepping

loglevel  = 1                       # amount of diagnostic output
				    
refine_grid = True                  # True to enable refinement

#################
m                  = gas.n_species

# specs index
ifuel = gas.species_index(fuel)
io2   = gas.species_index('O2')
in2   = gas.species_index('N2')

# Mixture composition
# Compute mixture composition based on equivalence ratio
gas.set_equivalence_ratio(phi, fuel, oxidizer)

x = np.zeros(m)

x[:] = gas.X  # This copies the mole fractions into x

################

#Set gas state to that of the unburned gas
gas.TPX = tin, p, x
    
#Create the free laminar premixed flame
f = FreeFlame(gas, initial_grid)
    
f.flame.set_steady_tolerances(default=tol_ss)
f.flame.set_transient_tolerances(default=tol_ts)
    
f.inlet.X = x
f.inlet.T = tin
    
f.transport_model = 'Mix'
    
#################################################################
# Program starts here
#################################################################
#First flame:
    
#No energy for starters
f.energy_enabled = False
    
#Refinement criteria
f.set_refine_criteria(ratio = 7.0, slope = 1, curve = 1)
    
#Max number of times the Jacobian will be used before it must be re-evaluated
f.set_max_jac_age(10, 10)
    
#Set time steps whenever Newton convergence fails
f.set_time_step(1.e-06, [1, 2, 5, 10]) #s
f.max_time_step_count = 3000
f.max_grid_points = 1000
    
#Calculation
f.solve(loglevel, refine_grid)
    
#################
#Second flame:
    
#Energy equation enabled
f.energy_enabled = True
    
#Refinement criteria when energy equation is enabled
f.set_refine_criteria(ratio = 5.0, slope = 0.5, curve = 0.5)
    
#Calculation 
f.solve(loglevel, refine_grid)
    
#################
#Third flame and so on ...:
f.set_refine_criteria(ratio = 3.0, slope = 0.1, curve = 0.1)    
f.solve(loglevel, refine_grid)
    
##################
##Fourth flame and so on ...:
f.set_refine_criteria(ratio = 2.0, slope = 0.05, curve = 0.05, prune = 0.01)    
f.solve(loglevel, refine_grid)
    
##################
##Fifth flame and so on ...
f.set_refine_criteria(ratio = 2.0, slope = 0.02, curve = 0.02, prune = 0.01)    
f.solve(loglevel, refine_grid)
    
print('mixture averaged flamespeed = ',f.velocity[0])



#Write the velocity, temperature, density, and mole fractions to a CSV file
# 6. Export results to CSV file
output_file = "onedim.csv"

nz = f.flame.n_points
with open(output_file, mode="w", newline="") as file:
    writer = csv.writer(file)
    # Write header
    header = ["X", "temp", "u", "rho"] + [f"X_{x}" for x in gas.species_names]
    writer.writerow(header)
    for kk in range(len(f.grid)):  
        f.set_gas_state(kk)    
        writer.writerow([f.grid[kk],gas.T, f.velocity[kk], gas.density]+list(gas.X))

print('export to CSV ...  DONE')



# PeleLM PMF file: MKS and mass fractions !
if (specformat == 'mole'):
    csv_file = str("pmf-X.txt")
    with open(csv_file, 'w') as outfile:
        writer = csv.writer(outfile, delimiter=' ', quotechar=' ', quoting=csv.QUOTE_MINIMAL)
        writer.writerow(['VARIABLES = "X" "temp" "u" "rho"']+['"X_'+x+'"' for x in gas.species_names])
        writer.writerow(['ZONE I=' + str(nz) + ' FORMAT=POINT' + ' SPECFORMAT=MOLE'])
        for kk in range(len(f.grid)):
            f.set_gas_state(kk)
            if (units == "MKS"): 
                writer.writerow([f.grid[kk],gas.T, f.velocity[kk], gas.density]+list(gas.X))
            else:
                writer.writerow([f.grid[kk]*1e2,gas.T, f.velocity[kk]*1e2, gas.density*1e-3]+list(gas.X))
if (specformat == 'mass'):
    csv_file = str("pmf-Y.txt")
    with open(csv_file, 'w') as outfile:
        writer = csv.writer(outfile, delimiter=' ', quotechar=' ', quoting=csv.QUOTE_MINIMAL)
        writer.writerow(['VARIABLES = "X" "temp" "u" "rho"']+['"Y_'+x+'"' for x in gas.species_names])
        writer.writerow(['ZONE I=' + str(nz) + ' FORMAT=POINT' + ' SPECFORMAT=MASS'])
        for kk in range(len(f.grid)):
            f.set_gas_state(kk)
            if (units == "MKS"): 
                writer.writerow([f.grid[kk],gas.T, f.velocity[kk], gas.density]+list(gas.Y))
            else:
                writer.writerow([f.grid[kk]*1e2,gas.T, f.velocity[kk]*1e2, gas.density*1e-3]+list(gas.Y))
