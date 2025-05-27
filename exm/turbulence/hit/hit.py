import numpy as np
import matplotlib.pyplot as plt
from numpy.fft import fftn, ifftn, fftfreq

def apply_band_limited_kolmogorov(k_mag, E_target=1.0, L_target=1.0):
    """Apply a band-limited Kolmogorov envelope centered around a length scale."""
    envelope = np.zeros_like(k_mag)

    k_peak = 1.0 / L_target
    k_min = k_peak / 2
    k_max = k_peak * 2

    k_min = 1
    k_max = 30

    print(" kmax kpeak kmin ",k_max,k_peak,k_min)


    mask = (k_mag >= k_min) & (k_mag <= k_max)
    #envelope[mask] = k_mag[mask]**(-5.0 / 6.0)

    envelope[mask] = k_mag[mask]**(-5.0 / 3.0)
    

    # Normalize to total target energy
    current_energy = np.sum(envelope**2)
    if current_energy > 0:
        envelope *= np.sqrt(E_target / current_energy)


    print(" current_energy ",current_energy)
    

    return envelope

def initialize_hit_velocity_field(N, L=1.0, L_t=1.0, seed=None, spectrum=None):
    if seed is not None:
        np.random.seed(seed)

    k = fftfreq(N, d=L/N) * 2 * np.pi
    kx, ky, kz = np.meshgrid(k, k, k, indexing='ij')
    k_squared = kx**2 + ky**2 + kz**2
    k_mag = np.sqrt(k_squared)
    k_squared[0, 0, 0] = 1.0  # avoid divide by zero

    # Random complex Fourier components
    fx = np.random.normal(size=(N, N, N)) + 1j * np.random.normal(size=(N, N, N))
    fy = np.random.normal(size=(N, N, N)) + 1j * np.random.normal(size=(N, N, N))
    fz = np.random.normal(size=(N, N, N)) + 1j * np.random.normal(size=(N, N, N))

    # Apply energy spectrum shaping if requested
    if spectrum == 'kolmogorov':
        envelope = np.zeros_like(k_mag)
        nonzero = k_mag > 0
        
        # Limit injection to a band of wavenumbers
        #envelope = apply_band_limited_kolmogorov(k_mag, E_target=1.0, L_target=L_t) 
    
        # Kolmogorov profile
        envelope[nonzero] = k_mag[nonzero]**(-5.0/3.0)
       
        fx *= envelope
        fy *= envelope
        fz *= envelope
    if spectrum == 'power_law':
        envelope = np.zeros_like(k_mag)
        nonzero = k_mag > 0     

        # all profile
        k0 = 2

        # Power (Garnier paper)  
        envelope[nonzero] = k_mag[nonzero]**4*np.exp(-2*(k_mag[nonzero]/k0)**2)
       
        fx *= envelope
        fy *= envelope
        fz *= envelope   

    # Enforce incompressibility (divergence-free)
    k_dot_f = kx * fx + ky * fy + kz * fz
    fx -= kx * k_dot_f / k_squared
    fy -= ky * k_dot_f / k_squared
    fz -= kz * k_dot_f / k_squared

    # Zero out mean mode
    fx[0, 0, 0] = fy[0, 0, 0] = fz[0, 0, 0] = 0.0

    # Transform to real space
    u = np.zeros((3, N, N, N), dtype=np.float64)
    u[0] = np.real(ifftn(fx))
    u[1] = np.real(ifftn(fy))
    u[2] = np.real(ifftn(fz))

    return u

def compute_energy_spectrum(u, L):
    N = u.shape[1]
    k = fftfreq(N, d=L/N) * 2 * np.pi
    kx, ky, kz = np.meshgrid(k, k, k, indexing='ij')
    k_mag = np.sqrt(kx**2 + ky**2 + kz**2)

    ux_hat = fftn(u[0])
    uy_hat = fftn(u[1])
    uz_hat = fftn(u[2])

    energy_density = 0.5 * (np.abs(ux_hat)**2 + np.abs(uy_hat)**2 + np.abs(uz_hat)**2)

    # Bin energy in shells of k
    k_max = np.max(k_mag)
    k_bins = np.linspace(0, k_max, N//2)
    spectrum = np.zeros_like(k_bins)

    shell_indices = np.digitize(k_mag.flatten(), k_bins)

    for i in range(1, len(k_bins)):
        shell_mask = (shell_indices == i)
        spectrum[i] = np.sum(energy_density.flatten()[shell_mask])

    k_centers = 0.5 * (k_bins[:-1] + k_bins[1:])
    return k_centers[1:], spectrum[1:-1]  # exclude k=0

def compute_integral_length_scale(k_vals, E_k):
    """Compute the integral length scale L_int from the energy spectrum E(k)."""
    dk = np.gradient(k_vals)
    numerator = np.sum(E_k / k_vals * dk)
    denominator = np.sum(E_k * dk)
    return numerator / denominator

def compute_taylor_microscale(k_vals, E_k):
    dk = np.gradient(k_vals)
    numerator = np.sum(E_k * dk)
    denominator = np.sum((k_vals**2) * E_k * dk)
    return np.sqrt(numerator / denominator)

def compute_rms_velocity(u):
    return np.sqrt(np.mean(u[0]**2 + u[1]**2 + u[2]**2))

def compute_taylor_reynolds_number(u, taylor_scale, nu):
    u_rms = compute_rms_velocity(u)
    return u_rms * taylor_scale / nu

def rescale_velocity_field_to_target_re_lambda(u, k_vals, E_k, nu, Re_lambda_target):
    # Step 1: Compute Taylor microscale
    lambda_taylor = compute_taylor_microscale(k_vals, E_k)
    
    # Step 2: Desired u_rms
    u_rms_target = Re_lambda_target * nu / lambda_taylor
    
    # Step 3: Current u_rms
    u_rms_actual = compute_rms_velocity(u)
    
    # Step 4: Rescale field
    scaling_factor = u_rms_target / u_rms_actual
    u_rescaled = u * scaling_factor
    
    return u_rescaled, u_rms_target, lambda_taylor

def compute_dissipation_rate(k_vals, E_k, nu):
    dk = np.gradient(k_vals)
    eps = 2 * nu * np.sum(k_vals**2 * E_k * dk)
    return eps

def compute_kolmogorov_length_scale(nu, eps):
    return (nu**3 / eps)**0.25

def export_velocity_field_binary(u, filename="velocity_field.bin"):
    """
    Export 3D velocity field (shape: [3, Nx, Ny, Nz]) to binary file.
    Header includes dimensions as 3 int32: Nx, Ny, Nz
    Data layout: ux, uy, uz in C-order (x-fastest).
    """
    Nx, Ny, Nz = u.shape[1:]
    with open(filename, "wb") as f:
        # Write header: grid size
        np.array([Nx, Ny, Nz], dtype=np.int32).tofile(f)
        # Write data: each component flattened in C-order
        for component in u:
            component.astype(np.float32).tofile(f)
    print(f"Exported velocity field to {filename}")


# --- Main ---
N = 64
L = 2 * np.pi
L_int0 = 0.3*L/10 #  Estimate of Integral Length scale
u_rms0 = 1
Re = 40   #  Initial estimate Reynolds based on Integral Length scale to compute viscosity
spectra_type = 'power_law'

nu = u_rms0*L_int0/Re

velocity_field = initialize_hit_velocity_field(N, L=L,L_t=L_int0, seed=42, spectrum=spectra_type)
k_vals, E_k = compute_energy_spectrum(velocity_field, L)

u_rms = compute_rms_velocity(velocity_field)
print(f"RMS velocity: u_rms bef scaling = {u_rms:.4e}")

# Rescale to match the prescribed rms0
velocity_field =  velocity_field*u_rms0/u_rms

#
u_rms = compute_rms_velocity(velocity_field)
print(f"RMS velocity: u_rms aft scaling = {u_rms:.4e}")

# Recompute spectrum after rescaling (optional)
k_vals, E_k = compute_energy_spectrum(velocity_field, L)

E_0 = E_k[1]

k_0 = 3**(-5/3) # for looking better in the plot

#print(f" E=0 ",E_0)


plt.figure(figsize=(6, 5))
plt.loglog(k_vals, E_k/E_0, label=' $E(k)$')
plt.loglog(k_vals, k_vals**(-5/3)/k_0, '--', label=r'$k^{-5/3}$ reference')
plt.xlabel('$k$')
plt.ylabel('$E(k)$')
plt.title('Normalized Energy Spectrum')
plt.legend()
plt.grid(True, which="both", ls="--")
plt.tight_layout()
plt.show()

# --- Compute and plot energy spectrum ---
dx = L/N
L_int = compute_integral_length_scale(k_vals, E_k)
lambda_taylor = compute_taylor_microscale(k_vals, E_k)
Re_lambda = compute_taylor_reynolds_number(velocity_field, lambda_taylor, nu)
Re_int = u_rms*L_int/nu
eps = compute_dissipation_rate(k_vals, E_k, nu)
eta = compute_kolmogorov_length_scale(nu, eps)

print("---------------------")
print(f"Re_lambda = {Re_lambda:.2f}, Re_int = {Re_int:.2f} , Re (defined) = {Re:.2f}")
print(f"viscosity: nu = {nu:.4e}")
print(f"Mean velocity: u_x = {np.mean(velocity_field[0]):.4e}")
print(f"RMS velocity: u_rms = {u_rms:.4e}")
print(f"Taylor-scale Reynolds number: Re_λ = {Re_lambda:.4f}")
print(f"Dissipation rate: ε = {eps:.4e}")
print(f"Integral length scale: L_int = {L_int:.4f}  L_int/L = {L_int/L:.4f} ")
print(f"Integral length scale: L_int/L_int0 =  {L_int/L_int0:.4f} ")
print(f"Taylor microscale: λ = {lambda_taylor:.4f}  λ/L = {lambda_taylor/L:.4f} ")
print(f"Kolmogorov length scale: η = {eta:.4e}  L/η = {L_int/eta:.4e}")
print(f"Max Resolution: dx = {dx:.4f}  L_int/dx = {L_int/dx:.4f} η/dx = {eta/dx:.4e}  ")

# PLOTS
# --- Plot slice of u_x ---
slice_index = N // 2
u_x_slice = velocity_field[0, :, :, slice_index]

plt.figure(figsize=(6, 5))
plt.imshow(u_x_slice, extent=[0, L, 0, L], origin='lower', cmap='viridis')
plt.colorbar(label='$u_x$')
plt.title('Slice of $u_x$ at $z = L/2$ with Kolmogorov Spectrum')
plt.xlabel('x')
plt.ylabel('y')
plt.tight_layout()
plt.show()

print(f" Exporting data ")
export_velocity_field_binary(velocity_field, filename="velocity_field.bin")

# plt.figure(figsize=(6, 5))
# plt.loglog(k_vals, E_k, label='Numerical $E(k)$')
# plt.loglog(k_vals, k_vals**(-5/3), '--', label=r'$k^{-5/3}$ reference')
# plt.xlabel('$k$')
# plt.ylabel('$E(k)$')
# plt.title('Energy Spectrum')
# plt.legend()
# plt.grid(True, which="both", ls="--")
# plt.tight_layout()
# plt.show()
