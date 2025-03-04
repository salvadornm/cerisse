import cantera as ct
import numpy as np
import matplotlib.pyplot as plt
import csv

# 1. Define gas mixture and mechanism
gas = ct.Solution("Li.yaml")  # Use a mechanism that includes H2/Air combustion

# 2. Set initial conditions
T0 = 1000.0  # Initial temperature [K]
P0 = ct.one_atm  # Initial pressure [Pa]
phi = 0.5  # Equivalence ratio (stoichiometric mixture)

# Define fuel and oxidizer mole fractions for a stoichiometric mixture
fuel = "H2"
oxidizer = {"O2": 1.0, "N2": 3.76}  # Air composition (O2/N2 ratio of 1:3.76)

# Compute mixture composition based on equivalence ratio
gas.set_equivalence_ratio(phi, fuel, oxidizer)
gas.TP = T0, P0  # Set temperature and pressure

# 3. Create a constant volume reactor
reactor = ct.IdealGasReactor(gas)               # constant volume reactor
#reactor = ct.IdealGasConstPressureReactor(gas) # constant pressure reactor
reactor_network = ct.ReactorNet([reactor])

# 4. Time integration setup
time = 0.0  # Start time [s]
end_time = 0.001  # End time [s]
times = []
temperatures = []

# 5. Integrate over time
while time < end_time:
    time = reactor_network.step()  # Integrate the system
    times.append(time)
    temperatures.append(reactor.T)


# 6. Export results to CSV file
output_file = "reactorV.csv"

with open(output_file, mode="w", newline="") as file:
    writer = csv.writer(file)
    writer.writerow(["Time (s)", "Temperature (K)"])  # Header
    for t, T in zip(times, temperatures):
        writer.writerow([t, T])

# 7. Plot results
plt.figure(figsize=(8, 5))
plt.plot(times, temperatures, label="Autoignition")
plt.xlabel("Time (s)")
plt.ylabel("Temperature (K)")
plt.title("H2/Air Autoignition in Cantera")
plt.legend()
plt.grid()
plt.show()

