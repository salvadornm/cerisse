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
x, t = sp.symbols('x t')
rho = 1 + 0.2 * sp.sin(2*sp.pi * x) 
u = 1
p = 1 + 0.3 * sp.cos(2*sp.pi * x)

# functions
rho_func = sp.lambdify(x, rho, modules=["numpy"])
u_func   = sp.lambdify(x, u, modules=["numpy"])
p_func   = sp.lambdify(x, p, modules=["numpy"])

x_vals = np.linspace(0, 1, 100)
rho_exact = rho_func(x_vals)
u_exact = u_func(x_vals)
p_exact = p_func(x_vals)

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

    xaxis = 0
    lineout = ds.ortho_ray(xaxis, (0, 0))
    srt = np.argsort(lineout["index", "x"])

    x_pred = np.array(lineout["index", "x"])[srt]
    rho_pred = np.array(lineout["boxlib", "Density"][srt])

    # Evaluate exact solution at numerical x points
    rho_exact_vals = rho_func(x_pred)

    # L2 norm of error
    abs_error = np.abs(rho_pred - rho_exact_vals)
    l2_error = np.sqrt(np.mean(abs_error**2))

    # Extract resolution (e.g., 16 from plot16)
    N = int(re.findall(r'\d+', dir_name)[0])
    resolutions.append(N)
    errors.append(l2_error)

    print(f"Resolution: {N}, L2 Error: {l2_error:.5e}")

# Step 3: Plot convergence
#resolutions = np.array(resolutions)
#errors = np.array(errors)

#plt.figure()
#plt.loglog(resolutions, errors, 'o-', label='L2 Error')
#plt.xlabel('Resolution (N)')
#plt.ylabel('L2 Error')
#plt.title('Grid Convergence Study')
#plt.grid(True, which='both')
#plt.legend()

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

