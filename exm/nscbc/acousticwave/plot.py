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


ds1 = yt.load("plot/plt00100")
ds2 = yt.load("plot/plt00200")
ds3 = yt.load("plot/plt00300")
ds4 = yt.load("plot/plt00400")


PREF=101325

print(' ##################################################### ')


xaxis = 0  # take a line cut along the x axis
lineout0 = ds0.ortho_ray(xaxis, (0, 0))   # cutting through the xaxis
lineout1 = ds1.ortho_ray(xaxis, (0, 0))   # cutting through the xaxis
lineout2 = ds2.ortho_ray(xaxis, (0, 0))   # cutting through the xaxis
lineout3 = ds3.ortho_ray(xaxis, (0, 0))   # cutting through the xaxis
lineout4 = ds4.ortho_ray(xaxis, (0, 0))   # cutting through the xaxis



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

srt2 = np.argsort(lineout2["index", "x"])
x2   = np.array(lineout2["index", "x"][srt2])
rho2 = np.array(lineout2["boxlib", "Density"][srt2])
p2   = np.array(lineout2["boxlib", "pressure"][srt2])
u2   = np.array(lineout2["boxlib", "x_velocity"][srt2])

srt3 = np.argsort(lineout3["index", "x"])
x3   = np.array(lineout3["index", "x"][srt3])
rho3 = np.array(lineout3["boxlib", "Density"][srt3])
p3   = np.array(lineout3["boxlib", "pressure"][srt3])
u3   = np.array(lineout3["boxlib", "x_velocity"][srt3])

srt4 = np.argsort(lineout4["index", "x"])
x4   = np.array(lineout4["index", "x"][srt4])
rho4 = np.array(lineout4["boxlib", "Density"][srt4])
p4   = np.array(lineout4["boxlib", "pressure"][srt4])
u4   = np.array(lineout4["boxlib", "x_velocity"][srt4])



# compute min/max 
rho1_min, rho1_max = np.min(rho1), np.max(rho1)
u1_min, u1_max = np.min(u1), np.max(u1)
u2_min, u2_max = np.min(u2), np.max(u2)
u3_min, u3_max = np.min(u3), np.max(u3)
u4_min, u4_max = np.min(u4), np.max(u4)

p1_min, p1_max = np.min(p1/PREF), np.max(p1/PREF)
p2_min, p2_max = np.min(p2/PREF), np.max(p2/PREF)
p3_min, p3_max = np.min(p3/PREF), np.max(p3/PREF)
p4_min, p4_max = np.min(p4/PREF), np.max(p4/PREF)


print(f" density: min={rho1_min:.4e}, max={rho1_max:.4e}")
print(f" pressure1: min={p1_min:.4e}, max={p1_max:.4e}")
print(f" pressure2: min={p2_min:.4e}, max={p2_max:.4e}")
print(f" pressure3: min={p3_min:.4e}, max={p3_max:.4e}")
print(f" pressure4: min={p4_min:.4e}, max={p4_max:.4e}")


print(f" velocity1: min={u1_min:.4e}, max={u1_max:.4e}")
print(f" velocity2: min={u2_min:.4e}, max={u2_max:.4e}")
print(f" velocity3: min={u3_min:.4e}, max={u3_max:.4e}")
print(f" velocity4: min={u4_min:.4e}, max={u4_max:.4e}")


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
axes[1].plot(x1, p1/PREF, label="t= 100", linestyle="--")
#axes[1].plot(x2, p2/PREF, label="t= 200", linestyle="--")
#axes[1].plot(x3, p3/PREF, label="t= 300", linestyle="--")
#axes[1].plot(x4, p4/PREF, label="t= 400", linestyle="-.")


axes[1].set_title("Pressure")
axes[1].set_xlabel("x")
axes[1].set_ylabel("p/p0")


axes[1].legend()


# Velocity plot
axes[2].plot(x0, u0, label="t=0")
axes[2].plot(x1, u1, label="t= 100", linestyle="--")
#axes[2].plot(x2, u2, label="t= 200", linestyle="--")
#axes[2].plot(x3, u3, label="t= 300", linestyle="--")
#axes[2].plot(x4, u4, label="t= 400", linestyle="-.")
axes[2].set_title("Velocity")
axes[2].set_xlabel("x")
axes[2].set_ylabel("Velocity")
axes[2].legend()

plt.tight_layout()
plt.show()

# Save the line plot into a file
#plot.save()


print(' .... DONE')
