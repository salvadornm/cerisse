import yt
import numpy as np
import matplotlib.pyplot as plt
from numpy.fft import fftn, fftfreq

def compute_energy_spectrum_from_yt(ds, vel_fields=("velocity_x", "velocity_y", "velocity_z")):
    ad = ds.all_data()
    
    # Get velocities as 3D numpy arrays
    ux = ad[vel_fields[0]].v  # .v gives NumPy array without units
    uy = ad[vel_fields[1]].v
    uz = ad[vel_fields[2]].v

    # Get the grid size
    N = ux.shape[0]
    assert ux.shape == uy.shape == uz.shape, "Velocity fields must have matching shapes"

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
    k_bins = np.linspace(0, k_max, N//2)
    spectrum = np.zeros(len(k_bins))

    shell_indices = np.digitize(k_mag, k_bins)

    for i in range(1, len(k_bins)):
        spectrum[i] = np.sum(energy_density_flat[shell_indices == i])

    # Midpoints of k bins for plotting
    k_centers = 0.5 * (k_bins[:-1] + k_bins[1:])
    
    return k_centers[1:], spectrum[1:-1]

# --- Usage ---
ds = yt.load("YourDataFile")  # or ds = yt.load(latest_file)
k_vals, E_k = compute_energy_spectrum_from_yt(ds)

# Plotting
plt.figure(figsize=(6, 5))
plt.loglog(k_vals, E_k, label='E(k)')
plt.loglog(k_vals, k_vals**(-5/3), '--', label=r'$k^{-5/3}$ ref')
plt.xlabel("Wavenumber k")
plt.ylabel("Energy Spectrum E(k)")
plt.title("3D Energy Spectrum from YT")
plt.grid(True, which="both", ls="--")
plt.legend()
plt.tight_layout()
plt.show()

