import numpy as np
import matplotlib.pyplot as plt
import glob
import pandas as pd

## load ref.data obtained from Cantera
df = pd.read_csv("reactorV.csv")
print(df.head())  # Show first 5 rows
# Extract columns
time = df["Time (s)"]
temperature = df["Temperature (K)"]
nsampling = 10
time_downsampled = df["Time (s)"][::nsampling]
temperature_downsampled = df["Temperature (K)"][::nsampling]


# ## load time.log obtained from Cerisse
ceris = pd.read_csv("probe.csv")
print(ceris.head())  # Show first 5 rows
# Extract column
time2 = ceris["time"]
temperature2 = ceris[" temperature((0)(15))"]

# Extract every 10th data point (nsampling)
nsampling = 50
time_downsampled2 = ceris["time"][::nsampling]
temperature_downsampled2 = ceris[" temperature((0)(15))"][::nsampling]

#Plot ref temperature vs time
plt.figure(figsize=(8, 5))

plt.plot(time_downsampled, temperature_downsampled, 'ko', label="Cantera") #every 10th 
plt.plot(time_downsampled2, temperature_downsampled2, label="Cerisse",color="r")  #every 50 th

plt.legend()

# Labels and title
plt.xlabel("Time (s)")
plt.ylabel("Temperature (K)")
plt.tight_layout()

plt.savefig("autoignition.png", dpi=300, bbox_inches='tight')  # Save with high resolution

plt.show()

#print(" ... Saved plot to autoignition!")
