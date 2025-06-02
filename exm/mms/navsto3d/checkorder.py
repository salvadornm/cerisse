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
gamma = sp.Rational(7, 5)  # gamma = 1.4

molweight = sp.Rational(29,1000)
Rgas =  8.31446261815324 / molweight
Cp   = gamma*Rgas / (gamma-1)

L_oo = 1
Re  = 1
Ma = 0.1
Pr = 0.72

rho_oo = 1.16
u_oo   = 152
v_oo   = 100
w_oo   = 0

csound = u_oo/Ma
p_oo   = rho_oo*csound*csound/gamma
T_oo   = p_oo/(Rgas*rho_oo)

visc = rho_oo*u_oo*L_oo/Re

cond = visc*Cp/Pr

# Exact solution
x, y, z, t = sp.symbols('x y z t')

## same as generate force
rho = rho_oo + 0.1 * sp.sin(2*sp.pi * x)  + 0.15 * sp.cos(2*sp.pi*y) + 0.2 * sp.sin(6*sp.pi*z) 
u = u_oo+ 27.0 * sp.sin(4*sp.pi * x) - 17.0 * sp.cos(2*sp.pi*y) +0.0 * sp.sin(4*sp.pi*z) 
v = v_oo + 69 *sp.sin(4*sp.pi*x)   + 0.0 * sp.cos(4*sp.pi*y)  + +0.0 * sp.sin(2*sp.pi*z) 
w = w_oo *sp.sin(4*sp.pi*x)   +0.0 * sp.cos(4*sp.pi*y)  + +0.0 * sp.sin(2*sp.pi*z)  
p = p_oo - 350 * sp.sin(2*sp.pi * x) + 60 * sp.cos(4*sp.pi *y) + 25 * sp.sin(6*sp.pi*z) 

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
errors_u =[]
errors_p =[]


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

    # Get predicted rho,u,p
    rho_pred = ad["boxlib", "Density"].ndarray_view()
    u_pred = ad["boxlib", "x_velocity"].ndarray_view()
    p_pred = ad["boxlib", "pressure"].ndarray_view()


    # Evaluate exact solution at those points
    rho_exact_pred = rho_func(X0, Y0, Z0)
    u_exact_pred   = u_func(X0, Y0, Z0)
    p_exact_pred   = p_func(X0, Y0, Z0)
    
    # Compute error density
    abs_error = np.abs(rho_pred - rho_exact_pred)
    l2_error = np.sqrt(np.mean(abs_error**2))
    errors.append(l2_error)
    
    abs_error_u = np.abs(u_pred - u_exact_pred)
    l2_error_u  = np.sqrt(np.mean(abs_error_u**2))
    errors_u.append(l2_error_u)

    abs_error_p = np.abs(p_pred - p_exact_pred)
    l2_error_p  = np.sqrt(np.mean(abs_error_p**2))
    errors_p.append(l2_error_p)


    # Extract resolution (e.g., 16 from plot16)
    N = int(re.findall(r'\d+', dir_name)[0])
    resolutions.append(N)

    



    print(f"Resolution: {N}, L2 Error: {l2_error:.5e}")


# Step 3: Plot convergence
resolutions = np.array(resolutions)
errors = np.array(errors)
errors_u = np.array(errors_u)
errors_p = np.array(errors_p)

e0   = errors[0]
e0_u = errors_u[0]
e0_p = errors_p[0]



plt.figure()
plt.loglog(resolutions, errors/e0, 'o-', label='rho')
plt.loglog(resolutions, errors_u/e0_u, 'o-', label='u')
plt.loglog(resolutions, errors_p/e0_p, 'o-', label='p')

# Add ideal convergence lines
N_ref = resolutions[0]
ideal_2nd = (resolutions / N_ref) ** -2
ideal_4th = (resolutions / N_ref) ** -4

plt.loglog(resolutions, ideal_2nd, 'k--', label='2nd-order')
plt.loglog(resolutions, ideal_4th, 'k-.', label='4th-order')

plt.xlabel('Resolution (N)')
plt.ylabel('L2 Error')
plt.title('Grid Convergence Study')
plt.grid(True, which='both')
plt.legend()

# Optional: fit and plot convergence rate
p = np.polyfit(np.log(resolutions), np.log(errors), 1)
order = -p[0]
#plt.text(resolutions[0], errors[0], f"Order ≈ {order:.2f}", fontsize=12)

print(' ##################################################### ')

print(f"Spatial Convergence Order : {order:.3f}")

plt.show()

