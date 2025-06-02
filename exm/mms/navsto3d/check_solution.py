import yt
import numpy as np
import matplotlib.pyplot as plt
import glob
import os
import sympy as sp
import re
import argparse




yt.funcs.mylog.setLevel("ERROR")  # or "CRITICAL" to suppress almost everything

#-------------------------------------------
# file to open data and plot 1D results
#-----------------------------------------


# Exact solution
x, t = sp.symbols('x t')
rho = 1 + 0.2 * sp.sin(2*sp.pi * x) 
u = sp.S(1)
p = 1 + 0.3 * sp.cos(2*sp.pi * x)

# functions
rho_func = sp.lambdify(x, rho, modules=["numpy"])
u_func   = sp.lambdify(x, u, modules=["numpy"])
p_func   = sp.lambdify(x, p, modules=["numpy"])

x_vals = np.linspace(0, 1, 100)
rho_exact = rho_func(x_vals)
u_exact = u_func(x_vals)
p_exact = p_func(x_vals)

parser = argparse.ArgumentParser(description="Check solution at a specific resolution.")
parser.add_argument("dir_name", type=str, help="Directory name (e.g., 'plot16')")
args = parser.parse_args()

dir_name = args.dir_name

file_name = f"{dir_name}/plt00001"

print(" plotting file .. \n",file_name)


datasets = [file_name]
ds = yt.load(datasets[0]);

xaxis = 0
lineout = ds.ortho_ray(xaxis, (0, 0))
srt = np.argsort(lineout["index", "x"])

x_pred   = np.array(lineout["index", "x"])[srt]
rho_pred = np.array(lineout["boxlib", "Density"][srt])
u_pred   = np.array(lineout["boxlib", "x_velocity"][srt])
p_pred   = np.array(lineout["boxlib", "pressure"][srt])


# Evaluate exact solution at numerical x points
rho_exact_vals = rho_func(x_pred)
u_exact_vals   = u_func(x_pred)
u_exact_vals = np.full_like(x_pred, u_exact_vals)
p_exact_vals   = p_func(x_pred)


# L2 norm of error
abs_error = np.abs(rho_pred - rho_exact_vals)
l2_error_R = np.sqrt(np.mean(abs_error**2))
print(f"L2 Density Error: {l2_error_R:.5e}")
abs_error = np.abs(u_pred - u_exact_vals)
l2_error_U = np.sqrt(np.mean(abs_error**2))
print(f"L2 Velocity Error: {l2_error_U:.5e}")
abs_error = np.abs(p_pred - p_exact_vals)
l2_error_P = np.sqrt(np.mean(abs_error**2))
print(f"L2 Pressure Error: {l2_error_P:.5e}")

# Plot predicted vs exact solution
# plt.figure()
# plt.plot(x_pred, rho_pred, '-', label='Numerical')        
# plt.plot(x_pred, rho_exact_vals, 'o', label='Exact')      
# plt.xlabel('x')
# plt.ylabel('Density')
# plt.title(f'Density Error: L2 = {l2_error_R:.5e}')
# plt.grid(True, which='both')
# plt.legend()
# plt.show()


fig, axs = plt.subplots(3, 1, figsize=(8, 10), sharex=True)

# # Plot 1: Density
axs[0].plot(x_pred, rho_pred, '-', label='Numerical')
axs[0].plot(x_pred, rho_exact_vals, 'o', label='Exact')
axs[0].set_ylabel('Density')
#axs[0].set_title('Density: Numerical vs Exact')
axs[0].set_title(f'Density Error: L2 = {l2_error_R:.5e}')
axs[0].grid(True)
axs[0].legend()

# # Plot 2: Velocity
axs[1].plot(x_pred, u_pred, '-', label='Numerical')
axs[1].plot(x_pred, u_exact_vals, 'o', label='Exact')
axs[1].set_ylabel('Velocity')
#axs[1].set_title('Velocity: Numerical vs Exact')
axs[2].set_title(f'Velocity Error: L2 = {l2_error_U:.5e}')
axs[1].grid(True)
#axs[1].legend()

# # Plot 3: Pressure
axs[2].plot(x_pred, p_pred, '-', label='Numerical')
axs[2].plot(x_pred, p_exact_vals, 'o', label='Exact')
axs[2].set_ylabel('Pressure')
axs[2].set_xlabel('x')
#axs[2].set_title('Pressure: Numerical vs Exact')
axs[2].set_title(f'Pressure Error: L2 = {l2_error_P:.5e}')
axs[2].grid(True)
#axs[2].legend()

plt.tight_layout()
plt.show()





print(' ##################################################### ')


plt.show()

