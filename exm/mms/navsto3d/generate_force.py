import sympy as sp
import numpy as np
from sympy import cxxcode
from datetime import datetime

# Symbols
x, y, z, t = sp.symbols('x y z t')
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

print(" csound = ", csound)
print(" visc = ", visc, " cond= ", cond)

print(" p_oo = ", p_oo, " T_oo = ",T_oo)


print(" Mach      = ", u_oo/csound)
print(" Reynolds  = ", rho_oo*u_oo*L_oo/visc)



# Manufactured solution
# steady
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

x0 = 0.0
y0 = 0.0
z0 = 0.0

# Define the domain
x_vals = np.linspace(0, 1, 100)
y_vals = np.linspace(0, 1, 100)
z_vals = np.linspace(0, 1, 100)

# Create 3D grid
X, Y, Z = np.meshgrid(x_vals, y_vals, z_vals, indexing='ij')

# Evaluate rho_func on the grid
rho_vals = rho_func(X, Y, Z)
u_vals   = u_func(X, Y, Z)
v_vals   = v_func(X, Y, Z)
w_vals   = w_func(X, Y, Z)
p_vals   = p_func(X, Y, Z)

# Compute max
rho_max  = np.max(rho_vals)
u_max    = np.max(u_vals)
v_max    = np.max(v_vals)
w_max    = np.max(w_vals)
p_max    = np.max(p_vals)

# print("Maximum of rho_func:", rho_max)
# print("Maximum of u_func:", u_max)
# print("Maximum of v_func:", v_max)
# print("Maximum of w_func:", w_max)
# print("Maximum of p_func:", p_max)
print("----------------------- \n")


gamma = 1.4

rho0 = rho_func(x0,y0,z0)
u0 = u_func(x0,y0,z0) 
v0 = v_func(x0,y0,z0)
w0 = w_func(x0,y0,z0)
p0 = p_func(x0,y0,z0)

vel =  np.sqrt( u0*u0 + v0*v0 +w0*w0 )

csound = np.sqrt(gamma * p0 / rho0)

print( " u0 v0 w0 ",u0,v0,w0)
print( " rho0 p0 ",rho0,p0)
print( " ufluc vfluc wfluc ",(u_max-u0)/u0,(v_max-v0)/v0,(w_max-w0)/(w0 + 0.000001) )
print( " rhofluc pfluc ",(rho_max-rho0)/rho0,(p_max-p0)/p0)
print( " csound=",csound, " vel=",vel)
print( " Mach= ",vel/csound)
print("----------------------- \n")

# Derived quantities
rho_u = rho * u
rho_v = rho * v
rho_w = rho * w
E = p / (gamma - 1) + 0.5 * rho * ( u**2 + v**2 + w**2)

T = p/(rho*Rgas)


# velocity gradients
dudx = sp.diff(u, x)
dudy = sp.diff(u, y)
dudz = sp.diff(u, z)
dvdx = sp.diff(v, x)
dvdy = sp.diff(v, y)
dvdz = sp.diff(v, z)
dwdx = sp.diff(w, x)
dwdy = sp.diff(w, y)
dwdz = sp.diff(w, z)

divu = sp.diff(u, x) + sp.diff(v, x) + sp.diff(w, z)

# shear stress
tau_xx = visc*(2*dudx - 2/3*divu) 
tau_xy = visc*(dudy + dvdx) 
tau_xz = visc*(dudz + dwdx)
tau_yy = visc*(2*dvdy - 2/3*divu) 
tau_yz = visc*(dvdz + dwdy)
tau_zz = visc*(2*dwdz - 2/3*divu) 

# Euler Fluxes x
F1x = rho_u
F2x = rho * u**2 + p
F3x = rho * u*v 
F4x = rho * u*w
F5x = u * (E + p)

# Visc Fluxes x
F1viscx = 0
F2viscx = tau_xx
F3viscx = tau_xy
F4viscx = tau_xz
F5viscx = cond*sp.diff(T, x) + u*tau_xx + v*tau_xy + w*tau_xz

# Euler Fluxes y
F1y = rho_v
F2y = rho * v*u
F3y = rho * v**2 + p 
F4y = rho * v*w
F5y = v * (E + p)

# Visc Fluxes y
F1viscy = 0
F2viscy = tau_xy
F3viscy = tau_yy
F4viscy = tau_yz
F5viscy = cond*sp.diff(T, y) + u*tau_xy + v*tau_yy + w*tau_yz

# Euler Fluxes z
F1z = rho_w
F2z = rho * w*u
F3z = rho * w*v
F4z = rho * w**2 + p
F5z = w * (E + p)

# Visc Fluxes y
F1viscz = 0
F2viscz = tau_xz
F3viscz = tau_yz
F4viscz = tau_zz
F5viscz = cond*sp.diff(T, z) + u*tau_xz + v*tau_yz + w*tau_zz

# Total fuxes x
Ft1x = F1x - F1viscx
Ft2x = F2x - F2viscx
Ft3x = F3x - F3viscx
Ft4x = F4x - F4viscx
Ft5x = F5x - F5viscx
# Total fuxes y
Ft1y = F1y - F1viscy
Ft2y = F2y - F2viscy
Ft3y = F3y - F3viscy
Ft4y = F4y - F4viscy
Ft5y = F5y - F5viscy
# Total fuxes z
Ft1z = F1z - F1viscz
Ft2z = F2z - F2viscz
Ft3z = F3z - F3viscz
Ft4z = F4z - F4viscz
Ft5z = F5z - F5viscz


# Spatial derivatives of fluxes
dF1 = sp.diff(Ft1x, x) + sp.diff(Ft1y, y) +  sp.diff(Ft1z, z) 
dF2 = sp.diff(Ft2x, x) + sp.diff(Ft2y, y) +  sp.diff(Ft2z, z) 
dF3 = sp.diff(Ft3x, x) + sp.diff(Ft3y, y) +  sp.diff(Ft3z, z) 
dF4 = sp.diff(Ft4x, x) + sp.diff(Ft4y, y) +  sp.diff(Ft4z, z) 
dF5 = sp.diff(Ft5x, x) + sp.diff(Ft5y, y) +  sp.diff(Ft5z, z) 

# Forcing terms: S = ∂U/∂t + ∂F/∂x -∂Fv/∂x
S1 = dF1
S2 = dF2
S3 = dF3
S4 = dF4
S5 = dF5

# Simplify
#S1_simplified = sp.simplify(S1)
#S2_simplified = sp.simplify(S2)
#S3_simplified = sp.simplify(S3)
#S4_simplified = sp.simplify(S4)
#S5_simplified = sp.simplify(S5)


# Output
print("rho :")
print(sp.latex(rho))

print("u :")
print(sp.latex(u))

print("v :")
print(sp.latex(v))

print("w :")
print(sp.latex(w))

print("p :")
print(sp.latex(p))


print("-----------------------\n")


# Output forcing terms
print("S_rho (mass equation: \n")
#sp.pprint(S1_simplified)
print(sp.latex(S1))

print("Sx_momentum (momentum equation): \n")
#sp.pprint(S2_simplified)
print(sp.latex(S2))

print("Sy_momentum (momentum equation): \n")
#sp.pprint(S3_simplified)
print(sp.latex(S3))

print("Sz_momentum (momentum equation): \n")
#sp.pprint(S4_simplified)
print(sp.latex(S4))

print("S_energy (energy equation):  \n")
#sp.pprint(S5_simplified)
print(sp.latex(S5))


print("-----------------------\n")

print(" Export to C++ \n")

# Write C++ header
hpp_content = f"""// This file is automatically generated
// Last update: {datetime.now().strftime("%d/%m/%Y %H:%M:%S")}
#ifndef MMS_HPP
#define MMS_HPP

#include <cmath>
#include <AMReX_REAL.H>

// Exact MMS solution //
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void mms_exact(const amrex::Real x, const amrex::Real y, const amrex::Real z, 
         amrex::Real& rho, amrex::Real& u, amrex::Real& v, amrex::Real& w, amrex::Real& p) {{
rho = {cxxcode(rho, standard='C++11')};
u = {cxxcode(u, standard='C++11')};
v = {cxxcode(v, standard='C++11')};
w = {cxxcode(w, standard='C++11')};
p = {cxxcode(p, standard='C++11')};
}}

// Source terms for the Euler equations //
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void mms_source(const amrex::Real x, const amrex::Real y, const amrex::Real z, 
         amrex::Real& Srho, amrex::Real& Srhou, amrex::Real& Srhov, amrex::Real& Srhow, amrex::Real& Srhoe) {{
Srho  = {cxxcode(S1, standard='C++11')};
Srhou = {cxxcode(S2, standard='C++11')};
Srhov = {cxxcode(S3, standard='C++11')};
Srhow = {cxxcode(S4, standard='C++11')};
Srhoe = {cxxcode(S5, standard='C++11')};
}}


#endif"""

with open("mms.h", "w") as file:
  file.write(hpp_content)

print("MMS exact solution and source term written to mms.h")


