import yt
import numpy as np
import matplotlib.pyplot as plt
import glob
import os
#-------------------------------------------
# file to open data and plot 1D results
# tst/tst1
#-----------------------------------------

print('One-dimensional advection of a density perturbation')
print(' opening data ..')

# load data
ds0 = yt.load("plot/plt00000")

list_of_files = glob.glob('./plot/*')
latest_file = max(list_of_files, key=os.path.getctime)
print(" LAST FILE=",latest_file)
ds = yt.load(latest_file)

PREF=101325

print(' ##################################################### ')


xaxis = 0  # take a line cut along the x axis
lineout0 = ds0.ortho_ray(xaxis, (0, 0))   # cutting through the xaxis
lineout1 = ds.ortho_ray(xaxis, (0, 0))   # cutting through the xaxis


# sort lineout0
srt = np.argsort(lineout0["index", "x"])
x0   = np.array(lineout0["index", "x"][srt])
rho0 = np.array(lineout0["boxlib", "Density"][srt])
p0   = np.array(lineout0["boxlib", "pressure"][srt])
u0   = np.array(lineout0["boxlib", "x_velocity"][srt])


# sort lineout1
srt1 = np.argsort(lineout1["index", "x"])
x1   = np.array(lineout1["index", "x"][srt1])
rho1 = np.array(lineout1["boxlib", "Density"][srt1])
p1   = np.array(lineout1["boxlib", "pressure"][srt1])
u1   = np.array(lineout1["boxlib", "x_velocity"][srt1])

# compute min/max 
rho1_min, rho1_max = np.min(rho1), np.max(rho1)
p1_min, p1_max = np.min(p1/PREF), np.max(p1/PREF)
u1_min, u1_max = np.min(u1), np.max(u1)

print(f" density: min={rho1_min:.4e}, max={rho1_max:.4e}")
print(f" pressure: min={p1_min:.4e}, max={p1_max:.4e}")
print(f" velocity: min={u1_min:.4e}, max={u1_max:.4e}")


# make plots (3 columns)
fig, axes = plt.subplots(1, 3, figsize=(16, 5))

# Density plot
axes[0].plot(x0, rho0, label="t=0")
axes[0].plot(x1, rho1, label="t=0.34", linestyle="--")
axes[0].set_title("Density")
axes[0].set_xlabel("x")
axes[0].set_ylabel("rho/rho0")
#axes[0].legend()

# Pressure plot
axes[1].plot(x0, p0/PREF, label="t=0")
axes[1].plot(x1, p1/PREF, label="t=0.34", linestyle="--")
axes[1].set_title("Pressure")
axes[1].set_xlabel("x")
axes[1].set_ylabel("p/p0")
#axes[1].legend()


# Velocity plot
axes[2].plot(x0, u0, label="t=0")
axes[2].plot(x1, u1, label="t=0.34", linestyle="--")
axes[2].set_title("Velocity")
axes[2].set_xlabel("x")
axes[2].set_ylabel("Velocity")
axes[2].legend()

plt.tight_layout()
plt.show()

# Save the line plot into a file
#plot.save()


print(' .... DONE')
