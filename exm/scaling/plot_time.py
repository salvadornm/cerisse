import numpy as np
import matplotlib.pyplot as plt
import argparse
import sys

# Parse command-line argument
parser = argparse.ArgumentParser(description="Plot a single quantity from time_probe.log")
parser.add_argument("field", type=str, help="Field to plot: kinetic, enstrophy, density, pressure")
args = parser.parse_args()
field = args.field.lower()

# Load data from CSV file
filename = "time_probe.log"
try:
    data = np.loadtxt(filename, delimiter=',', skiprows=1)
except Exception as e:
    print(f"Error reading {filename}: {e}")
    exit(1)


# Map field names to column index and labels
field_map = {
    "kinetic": (1, "Kinetic Energy"),
    "enstrophy": (2, "Enstrophy"),
    "density": (3, "Density"),
    "pressure": (4, "Pressure")
}

# Unpack columns
if field not in field_map:
    print(f"Unknown field '{field}'. Choose from: {', '.join(field_map.keys())}")
    sys.exit(1)

col_index, label = field_map[field]

# Extract time and selected field
time = data[:, 0]
quantity = data[:, col_index]

# non-dimensional time
t0 = 0.02879768262 # from screen L/Urms

# normalised
q0 = quantity[0]


# Plot
plt.figure(figsize=(8, 5))
plt.plot(time/t0, quantity/q0, label=label)
plt.xlabel("Time")
plt.ylabel(label)
plt.title(f"{label} vs Time")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()
