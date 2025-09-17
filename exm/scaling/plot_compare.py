import numpy as np
import matplotlib.pyplot as plt
import argparse



file1 = "refdata/cerisse_teno6_128.log"

col_probe = 2 # enstrophy column

t0 = 0.02879768262 # from screen L/Urms

###########
file2 = "refdata/spectral_Re1600_512.gdiag"
col_spectral= 3 # enstrophy column

####
file3 = "refdata/opensbli_128.dat"
col_sbli= 2

## current data
file0 = "time_probe.log"

# --- Load file0 (timeprob.log) ---
data0 = np.loadtxt(file0, delimiter=',', skiprows=1)
time0 = data0[:, 0]
enstrophy0 = data0[:, col_probe]

# --- Load file1 (time.log) ---
data1 = np.loadtxt(file1, delimiter=',', skiprows=1)
time1 = data1[:, 0]
enstrophy1 = data1[:, col_probe]

# --- Load file2 (spectral.gdiag) ---
# Expected format: Time  Energy  -dE/dt  Enstrophy
data2 = np.loadtxt(file2)
time2 = data2[:, 0]
enstrophy2 = data2[:, col_spectral]

# --- Load file3 (sbli.dat) ---
# Expected format: Time  Energy  Enstrophy
data3 = np.loadtxt(file3)
time3 = data3[:, 0]
enstrophy3 = data3[:, col_sbli]


E00 = enstrophy0[0]
E10 = enstrophy1[0]
E20 = enstrophy2[0]
E30 = enstrophy3[0]

plt.figure(figsize=(9, 5))


# Plot cerisse ref ---
plt.plot(time1/t0, enstrophy1/E10, label="Cerisse TENO6", linewidth=2)

# Plot spectral: faint line + markers every 20 points
plt.plot(time2, enstrophy2 / E20, color='gray', linewidth=1, label='_nolegend_')
plt.plot(time2[::20], (enstrophy2 / E20)[::20], 'o', color='black', markersize=4, label="Reference")

# Plot OpenSBLI
plt.plot(time3, enstrophy3/E30, label="OpenSBLI", linestyle='-', linewidth=2)

# Plot current data
plt.plot(time0/t0, enstrophy0/E00, label="current", linestyle='-', linewidth=2)


plt.xlim(0, 20)

plt.xlabel("Time")
plt.ylabel("Enstrophy")
plt.title("Normalised Enstrophy ")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

