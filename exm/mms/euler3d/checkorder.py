import yt
import numpy as np
import matplotlib.pyplot as plt
import glob
import os
import sympy as sp
import re



yt.funcs.mylog.setLevel("ERROR")  # or "CRITICAL" to suppress almost everything

#-------------------------------------------
# file to open data and plot 1D results
#-----------------------------------------


# Exact solution
x, y, z, t = sp.symbols('x y z t')
rho = 1.16 + 0.1 * sp.sin(2*sp.pi * x)  + 0.15 * sp.cos(2*sp.pi*y) + 0.2 * sp.sin(6*sp.pi*z) 
u = 152+ 27.0 * sp.sin(4*sp.pi * x) - 17.0 * sp.cos(2*sp.pi*y) +0.0 * sp.sin(4*sp.pi*z) 
v = 100 + 69 *sp.sin(4*sp.pi*x)   + 0.0 * sp.cos(4*sp.pi*y)  + +0.0 * sp.sin(2*sp.pi*z) 
w = 0.0 *sp.sin(4*sp.pi*x)   +0.0 * sp.cos(4*sp.pi*y)  + +0.0 * sp.sin(2*sp.pi*z)  
p = 1e5 - 350 * sp.sin(2*sp.pi * x) + 60 * sp.cos(4*sp.pi *y) + 25 * sp.sin(6*sp.pi*z) 


# functions
rho_func = sp.lambdify((x, y, z), rho, modules=["numpy"])
u_func   = sp.lambdify((x, y, z),   u, modules=["numpy"])
v_func   = sp.lambdify((x, y, z),   v, modules=["numpy"])
w_func   = sp.lambdify((x, y, z),   w, modules=["numpy"])
p_func   = sp.lambdify((x, y, z),   p, modules=["numpy"])


# Define the domain
x_vals = np.linspace(0, 1, 100)
y_vals = np.linspace(0, 1, 100)
z_vals = np.linspace(0, 1, 100)

# Create 3D grid
X, Y, Z = np.meshgrid(x_vals, y_vals, z_vals, indexing='ij')

# Evaluate rho_func on the grid
rho_exact = rho_func(X, Y, Z)
u_exact   = u_func(X, Y, Z)
v_exact   = v_func(X, Y, Z)
w_exact   = w_func(X, Y, Z)
p_exact   = p_func(X, Y, Z)

# Step 2: Define directories
base_dirs = sorted([d for d in os.listdir() if re.match(r'plot\d+', d)],
                   key=lambda s: int(re.findall(r'\d+', s)[0]))

resolutions = []
errors = []


for dir_name in base_dirs:
    #print(f"Processing {dir_name} ...")

    file_name = f"{dir_name}/plt00001"

    datasets = [file_name]
    ds = yt.load(datasets[0]);

    # Get all cell centers
    ad = ds.all_data()
    X0 = ad["index", "x"].to("code_length").ndarray_view()
    Y0 = ad["index", "y"].to("code_length").ndarray_view()
    Z0 = ad["index", "z"].to("code_length").ndarray_view()

    # Get predicted density
    rho_pred = ad["boxlib", "Density"].ndarray_view()

    # Evaluate exact solution at those points
    rho_exact_pred = rho_func(X0, Y0, Z0)

    # Compute error density
    abs_error = np.abs(rho_pred - rho_exact_pred)
    l2_error = np.sqrt(np.mean(abs_error**2))

    # Extract resolution (e.g., 16 from plot16)
    N = int(re.findall(r'\d+', dir_name)[0])
    resolutions.append(N)
    errors.append(l2_error)

    print(f"Resolution: {N}, L2 Error: {l2_error:.5e}")

# Step 3: Plot convergence
resolutions = np.array(resolutions)
errors = np.array(errors)

plt.figure()
plt.loglog(resolutions, errors, 'o-', label='L2 Error')
plt.xlabel('Resolution (N)')
plt.ylabel('L2 Error')
plt.title('Grid Convergence Study')
plt.grid(True, which='both')
plt.legend()

# Step 3: Plot convergence
resolutions = np.array(resolutions)
errors = np.array(errors)

plt.figure()
plt.loglog(resolutions, errors, 'o-', label='L2 Error')
plt.xlabel('Resolution (N)')
plt.ylabel('L2 Error')
plt.title('Grid Convergence Study')
plt.grid(True, which='both')
plt.legend()

# Optional: fit and plot convergence rate
p = np.polyfit(np.log(resolutions), np.log(errors), 1)
order = -p[0]
plt.text(resolutions[0], errors[0], f"Order ≈ {order:.2f}", fontsize=12)

print(' ##################################################### ')

print(f"Spatial Convergence Order : {order:.3f}")

plt.show()

