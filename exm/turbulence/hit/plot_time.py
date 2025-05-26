import numpy as np
import matplotlib.pyplot as plt

# Load data from CSV file
filename = "time_probe.log"
try:
    data = np.loadtxt(filename, delimiter=',', skiprows=1)
except Exception as e:
    print(f"Error reading {filename}: {e}")
    exit(1)

# Unpack columns
time = data[:, 0]
kinetic_energy = data[:, 1]
enstrophy = data[:, 2]
density = data[:, 3]
pressure = data[:, 4]

# Plotting
plt.figure(figsize=(10, 6))
plt.plot(time, kinetic_energy, label="Kinetic Energy")
#plt.plot(time, enstrophy, label="Enstrophy")
#plt.plot(time, pressure, label="Pressure", linestyle='--')
#plt.plot(time, density, label="Density", linestyle=':')

plt.xlabel("Time")
plt.ylabel("Quantity")
plt.title("Time Evolution of Physical Quantities")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()

