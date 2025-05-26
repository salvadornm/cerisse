# python script to exctact spectra

import yt
import numpy as np
import matplotlib.pyplot as plt
from numpy.fft import fftn, fftfreq
import glob
import os
import argparse


def compute_energy_spectrum_from_yt(ds, vel_fields=("velocity_x", "velocity_y", "velocity_z")):
    ad = ds.all_data()
    
    # Get velocities as 3D numpy arrays
    ux = ad[vel_fields[0]].v  # .v gives NumPy array without units
    uy = ad[vel_fields[1]].v
    uz = ad[vel_fields[2]].v

    # Get the grid size
    #N = ux.shape[0]
    #assert ux.shape == uy.shape == uz.shape, "Velocity fields must have matching shapes"

    grid_shape = ds.domain_dimensions
    assert grid_shape[0] == grid_shape[1] == grid_shape[2], "Grid must be cubic"
    N = grid_shape[0]

    print(" mesh N = ",N)

    

    # FFT of each component
    ux_hat = fftn(ux)
    uy_hat = fftn(uy)
    uz_hat = fftn(uz)

    # Compute energy density in spectral space
    energy_density = 0.5 * (np.abs(ux_hat)**2 + np.abs(uy_hat)**2 + np.abs(uz_hat)**2)

    # Create wavenumber grid
    L = ds.domain_width[0].v  # Assuming cubic box
    k = fftfreq(N, d=L/N) * 2 * np.pi
    kx, ky, kz = np.meshgrid(k, k, k, indexing='ij')
    k_mag = np.sqrt(kx**2 + ky**2 + kz**2).flatten()

    # Bin energy into spherical shells
    energy_density_flat = energy_density.flatten()
    k_max = np.max(k)
    # linear
    k_bins = np.linspace(0, k_max, N//2)
    # log
    #k_bins = np.logspace(np.log10(np.min(k[k > 0])), np.log10(k_max), num=50)

    # Bin centers (length = num_bins - 1)
    k_centers = 0.5 * (k_bins[:-1] + k_bins[1:])

    # Spectrum computed over len(k_bins) bins
    spectrum = np.zeros(len(k_bins) - 1)

    # Assign energy to spectrum bins (matching length of k_centers)
    shell_indices = np.digitize(k_mag, k_bins) - 1  # subtract 1 to match spectrum index
    for i in range(len(spectrum)):
        spectrum[i] = np.sum(energy_density_flat[shell_indices == i])

    # Filter out any zero entries (optional, for clean plotting)
    valid = (k_centers > 0) & (spectrum > 0)

    # Return safely aligned arrays
    return k_centers[valid], spectrum[valid]

# --- Usage ---


# take dirname as argument
parser = argparse.ArgumentParser(description="Check solution at a specific resolution.")
parser.add_argument("dir_name", type=str, help="Directory name (e.g., 'plot16')")
args = parser.parse_args()

dir_name = args.dir_name

# list of files in dir
filelist_name = f"{dir_name}/*"
list_of_files = glob.glob(filelist_name)
latest_file = max(list_of_files, key=os.path.getctime)
print(" LAST FILE=",latest_file)

ds = yt.load(latest_file)
###############################

k_vals, E_k = compute_energy_spectrum_from_yt(ds)

E_0 = E_k[1]
k_0 = 1**(-5/3) # for Kolmogorov looking better in the plot
# Plotting
plt.figure(figsize=(6, 5))
# plt.loglog(k_vals, E_k/E_0, label='E(k)')
# plt.loglog(k_vals, k_vals**(-5/3)/k_0, '--', label=r'$k^{-5/3}$ ref')
plt.loglog(k_vals[1:], (E_k/E_0)[1:], label='E(k)')
plt.loglog(k_vals[1:], k_vals[1:]**(-5/3)/k_0, '--', label=r'$k^{-5/3}$ ref')


plt.xlabel("Wavenumber k")
plt.ylabel("Energy Spectrum E(k)")
plt.title("3D Energy Spectrum from YT")
plt.grid(True, which="both", ls="--")
plt.legend()
plt.tight_layout()
plt.show()

