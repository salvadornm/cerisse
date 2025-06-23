import sympy as sp
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import derivatives
import visco_solver


FU = sp.FU


NavierStokes = 0  # if 1 solve Navier Stokes
Model        = 2  # 0:central  1:conservative flux/FD  2:conservative flux/FV
FiniteVolume = 1  # 0:FD      1: FV
ErrorMeasure   = 2  # 1:L1  2:L2  otherse Loo
EulerFlux      = 1  # 1:fluxes are on, otherwise off
ErrorQ         = 1  # 0: error in U  1: error in Q
# Symbols
x, y, z, t = sp.symbols('x y z t')
gamma = sp.Rational(7, 5)  # gamma = 1.4

molweight = sp.Rational(29,1000)
Rgas =  8.31446261815324 / molweight
Cp   = gamma*Rgas / (gamma-1)
o_Rgas = 1.0/Rgas

L_oo = 1
Re  = 100 #
Ma = 0.1
Pr = 0.72

rho_oo = 1.16
u_oo   = 152
v_oo   = 100
w_oo   = 0.0

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
# rho = rho_oo + 0.1 * sp.sin(2*sp.pi * x)  + 0.15 * sp.cos(2*sp.pi*y) + 0.2 * sp.sin(6*sp.pi*z) 
# u = u_oo+ 27.0 * sp.sin(4*sp.pi * x) - 17.0 * sp.cos(2*sp.pi*y) +0.0 * sp.sin(4*sp.pi*z) 
# v = v_oo + 69 *sp.sin(4*sp.pi*x)   + 0.0 * sp.cos(4*sp.pi*y)  + +0.0 * sp.sin(2*sp.pi*z) 
# w = w_oo *sp.sin(4*sp.pi*x)   +0.0 * sp.cos(4*sp.pi*y)  + +0.0 * sp.sin(2*sp.pi*z)  
# p = p_oo - 350 * sp.sin(2*sp.pi * x) + 60 * sp.cos(4*sp.pi *y) + 25 * sp.sin(6*sp.pi*z) 

# snm
rho = rho_oo + 0.1 * rho_oo* sp.sin(2*sp.pi * x)
u   = u_oo +  0.1*u_oo*sp.sin(2*sp.pi * y) + 0.01*u_oo*sp.sin(2*sp.pi * z) 
v   = 0.01*u_oo *sp.sin(4*sp.pi*x)  
w   = 0.0
p   = p_oo  + 0.1*p_oo * sp.sin(2*sp.pi * x) 


rho = rho_oo
p   = p_oo


# functions
rho_func = sp.lambdify((x, y, z), rho, modules=["numpy"])
u_func   = sp.lambdify((x, y, z),   u, modules=["numpy"])
v_func   = sp.lambdify((x, y, z),   v, modules=["numpy"])
w_func   = sp.lambdify((x, y, z),   w, modules=["numpy"])
p_func   = sp.lambdify((x, y, z),   p, modules=["numpy"])

# Integrate function symbolic to get rhof uf vf wf pf
# FINTE VOLUME   ceate function that computes rhof = 1/(dx*dy*dz) I(x+dx,y+) 

dx, dy, dz = sp.symbols('dx dy dz', positive=True)
# Define integration bounds
x_bounds = (x, x - dx/2, x + dx/2)
y_bounds = (y, y - dy/2, y + dy/2)
z_bounds = (z, z - dz/2, z + dz/2)

# Compute the integral
if (FiniteVolume==1):
    rhoI = sp.integrate(rho, z_bounds, y_bounds, x_bounds)
    uI   = sp.integrate(u  , z_bounds, y_bounds, x_bounds)
    vI   = sp.integrate(v  , z_bounds, y_bounds, x_bounds)
    wI   = sp.integrate(w  , z_bounds, y_bounds, x_bounds)
    pI   = sp.integrate(p  , z_bounds, y_bounds, x_bounds)

    rhoI_func = sp.lambdify((x, y, z, dx, dy, dz), rhoI, modules=["numpy"])
    uI_func   = sp.lambdify((x, y, z, dx, dy, dz), uI, modules=["numpy"])
    vI_func   = sp.lambdify((x, y, z, dx, dy, dz), vI, modules=["numpy"])
    wI_func   = sp.lambdify((x, y, z, dx, dy, dz), wI, modules=["numpy"])
    pI_func   = sp.lambdify((x, y, z, dx, dy, dz), pI, modules=["numpy"])


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
rho_min  = np.min(rho_vals)
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

csound0 = np.sqrt(gamma * p0 / rho0)

# print( " u0 v0 w0 ",u0,v0,w0)
# print( " rho0 p0 ",rho0,p0)
# print( " ufluc vfluc wfluc ",(u_max-u0)/u0,(v_max-v0)/v0,(w_max-w0)/(w0 + 0.000001) )
# print( " rhofluc pfluc ",(rho_max-rho0)/rho0,(p_max-p0)/p0)
print( " csound csound0=",csound,csound0, " vel=",vel)
print( " Mach= ",u_oo/csound)
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

divu = dudx+ dvdy + dwdz

# shear stress
tau_xx = visc*(2*dudx - 2/3*divu) 
tau_xy = visc*(dudy + dvdx) 
tau_xz = visc*(dudz + dwdx)
tau_yy = visc*(2*dvdy - 2/3*divu) 
tau_yz = visc*(dvdz + dwdy)
tau_zz = visc*(2*dwdz - 2/3*divu) 

if (EulerFlux ==1):
    # Euler Fluxes x
    F1x = rho_u
    F2x = rho * u**2 + p
    F3x = rho * u*v 
    F4x = rho * u*w
    F5x = u * (E + p)
    # Euler Fluxes y
    F1y = rho_v
    F2y = rho * v*u
    F3y = rho * v**2 + p 
    F4y = rho * v*w
    F5y = v * (E + p)
    # Euler Fluxes z
    F1z = rho_w
    F2z = rho * w*u
    F3z = rho * w*v
    F4z = rho * w**2 + p
    F5z = w * (E + p)
else:
    F1x = 0
    F2x = 0
    F3x = 0
    F4x = 0
    F5x = 0
    F1y = 0
    F2y = 0
    F3y = 0
    F4y = 0
    F5y = 0
    F1z = 0
    F2z = 0
    F3z = 0
    F4z = 0
    F5z = 0 

if NavierStokes  >0:
    # Visc Fluxes x
    F1viscx = 0
    F2viscx = tau_xx
    F3viscx = tau_xy
    F4viscx = tau_xz
    F5viscx = cond*sp.diff(T, x) + u*tau_xx + v*tau_xy #+ w*tau_xz
    # Visc Fluxes y
    F1viscy = 0
    F2viscy = tau_xy
    F3viscy = tau_yy
    F4viscy = tau_yz
    F5viscy = cond*sp.diff(T, y) + u*tau_xy + v*tau_yy #+ w*tau_yz
    # Visc Fluxes z
    F1viscz = 0
    F2viscz = tau_xz
    F3viscz = tau_yz
    F4viscz = tau_zz
    F5viscz = cond*sp.diff(T, z) + u*tau_xz + v*tau_yz #+ w*tau_zz
else:
    F1viscx = 0
    F2viscx = 0
    F3viscx = 0
    F4viscx = 0
    F5viscx = 0
    F1viscy = 0
    F2viscy = 0
    F3viscy = 0
    F4viscy = 0
    F5viscy = 0
    F1viscz = 0
    F2viscz = 0
    F3viscz = 0
    F4viscz = 0
    F5viscz = 0
    
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

#print(sp.latex(S2))

#print(sp.latex(S5))

# functions
S1_func   = sp.lambdify((x, y, z), S1, modules=["numpy"])
S2_func   = sp.lambdify((x, y, z), S2, modules=["numpy"])
S3_func   = sp.lambdify((x, y, z), S3, modules=["numpy"])
S4_func   = sp.lambdify((x, y, z), S4, modules=["numpy"])
S5_func   = sp.lambdify((x, y, z), S5, modules=["numpy"])

# Compute the integral of the source term
if (FiniteVolume==1):
    print(" Computing integrals..")
    S1I   = sp.integrate(S1  , z_bounds, y_bounds, x_bounds)
    print(" S1 ... DONE")
    S2I   = sp.integrate(S2  , z_bounds, y_bounds, x_bounds)
    print(" S2 ... DONE")
    S3I   = sp.integrate(S3  , z_bounds, y_bounds, x_bounds)
    print(" S3 ... DONE")
    S4I   = sp.integrate(S4  , z_bounds, y_bounds, x_bounds)
    print(" S4 ... DONE")

   # S5aux = FU['TR8'](S5)

    # print(" FU... DONE")
    # print(sp.latex(S5aux ))

    S5I   = sp.integrate( S5, z_bounds, y_bounds, x_bounds)
    #S5I = S5 #temp snm

    #S5_simplified = sp.simplify(S5)
    #print(sp.latex(S5_simplified)
    
    print(" S5 ... DONE")


    S1I_func = sp.lambdify((x, y, z, dx, dy, dz), S1I, modules=["numpy"])
    S2I_func = sp.lambdify((x, y, z, dx, dy, dz), S2I, modules=["numpy"])
    S3I_func = sp.lambdify((x, y, z, dx, dy, dz), S3I, modules=["numpy"])
    S4I_func = sp.lambdify((x, y, z, dx, dy, dz), S4I, modules=["numpy"])
    S5I_func = sp.lambdify((x, y, z, dx, dy, dz), S5I, modules=["numpy"])

    # S1I_func = sp.lambdify((x, y, z, dx, dy, dz), S1, modules=["numpy"])
    # S2I_func = sp.lambdify((x, y, z, dx, dy, dz), S2, modules=["numpy"])
    # S3I_func = sp.lambdify((x, y, z, dx, dy, dz), S3, modules=["numpy"])
    # S4I_func = sp.lambdify((x, y, z, dx, dy, dz), S4, modules=["numpy"])
    # S5I_func = sp.lambdify((x, y, z, dx, dy, dz), S5, modules=["numpy"])





# print(" test ...")
# ds = 0.01
# x1 = 0.5
# y1 = 0.5
# z1 = 0.5
# dvol = ds*ds*ds

# S1_val = S1I.evalf(subs={x: x1, y: y1, z: z1, dx: ds, dy: ds, dz: ds})
# print("Srho FINITE VOLUME     Sf  = 1/V int S dV:", S1_val/dvol)
# print("Srho FINITE DIFFERENCE S   = S(x,y,z    ):", S1_func(x1,y1,z1))
# print("Srho FINITE VOLUME     Sf_func  :", S1I_func(x1,y1,z1,ds,ds,ds)/dvol)
# S2_val = S2I.evalf(subs={x: x1, y: y1, z: z1, dx: ds, dy: ds, dz: ds})
# print("Su FINITE VOLUME     Sf  = 1/V int S dV:", S2_val/dvol)
# print("Su FINITE DIFFERENCE S   = S(x,y,z    ):", S2_func(x1,y1,z1))
# print("Su FINITE VOLUME     Sf_func  :", S2I_func(x1,y1,z1,ds,ds,ds)/dvol)
# S3_val = S3I.evalf(subs={x: x1, y: y1, z: z1, dx: ds, dy: ds, dz: ds})
# print("Sv FINITE VOLUME     Sf  = 1/V int S dV:", S3_val/dvol)
# print("Sv FINITE DIFFERENCE S   = S(x,y,z    ):", S3_func(x1,y1,z1))
# print("Sv FINITE VOLUME     Sf_func  :", S3I_func(x1,y1,z1,ds,ds,ds)/dvol)
# S4_val = S4I.evalf(subs={x: x1, y: y1, z: z1, dx: ds, dy: ds, dz: ds})
# print("Sw FINITE VOLUME     Sf  = 1/V int S dV:", S4_val/dvol)
# print("Sw FINITE DIFFERENCE S   = S(x,y,z    ):", S4_func(x1,y1,z1))
# print("Sw FINITE VOLUME     Sf_func  :", S4I_func(x1,y1,z1,ds,ds,ds)/dvol)
# S5_val = S5I.evalf(subs={x: x1, y: y1, z: z1, dx: ds, dy: ds, dz: ds})
# print("Se FINITE VOLUME     Sf  = 1/V int S dV:", S5_val/dvol)
# print("Se FINITE DIFFERENCE S   = S(x,y,z    ):", S5_func(x1,y1,z1))
# print("Se FINITE VOLUME     Sf_func  :", S5I_func(x1,y1,z1,ds,ds,ds)/dvol)


#-------- create solver

Lx = L_oo
Ly = L_oo
Lz = L_oo


order = derivatives.order  #<============ 
HALO = int(order/2)


# Resolutions to test (finer will not work)
resolutions = [8, 16, 32, 64]
#resolutions = [8, 16, 32]
#resolutions = [4]

errors_rho  = []
errors_u    = []
errors_v    = []
errors_w    = []
errors_p    = []

NCOMP = 5


print(" .... START CALCULATING .... ")
############  loop over resolutions
for N in resolutions:

    Nx = N
    Ny = N
    Nz = N
    dxx = Lx/Nx
    dyy = Ly/Ny
    dzz = Lz/Nz

    o_dx = 1.0/dxx
    o_dy = 1.0/dyy
    o_dz = 1.0/dzz
    o_dvol = o_dx*o_dy*o_dz

    #  array   0,1, ... HALO -1       N + HALO -1
    #  64 cells  HALO 2    0,1, HALO...... 64,
    ist = order         # 2/4/6 (first inner point)
    ien = ist + Nx - 1  # 8/9/10
    jst = order        
    jen = ist + Ny - 1
    kst = order          
    ken = kst + Nz - 1

    index_start = [ist,jst,kst]
    index_end   = [ien,jen,ken]

    HBC = HALO
    index_start_bc = [ist-HBC,jst-HBC,kst-HBC]  #1/2/3
    index_end_bc   = [ien+HBC,jen+HBC,ken+HBC]

    NTX = Nx + 2*order  # originally 2*HALO 
    NTY = Ny + 2*order
    NTZ = Nz + 2*order

    # print( "ist jst kst ", ist, jst, kst)
    # print( "ien jen ken ", ien, jen, ken)
    # print( "NTX NTY NTZ ", NTX, NTY, NTZ)

    # 3D array with all elements initialized to 0.0
    RHOU = np.zeros((NTX, NTY, NTZ))  
    RHOV = np.zeros((NTX, NTY, NTZ))  
    RHOW = np.zeros((NTX, NTY, NTZ))  
    RHOET = np.zeros((NTX, NTY, NTZ))  

    RHOU0 = np.zeros((NTX, NTY, NTZ))  
    RHOV0 = np.zeros((NTX, NTY, NTZ))  
    RHOW0 = np.zeros((NTX, NTY, NTZ))  
    RHOET0 = np.zeros((NTX, NTY, NTZ))  

    # errors
    ERR_RHO  = np.zeros((NTX, NTY, NTZ))  
    ERR_U    = np.zeros((NTX, NTY, NTZ))  
    ERR_V    = np.zeros((NTX, NTY, NTZ))  
    ERR_W    = np.zeros((NTX, NTY, NTZ))  
    ERR_P    = np.zeros((NTX, NTY, NTZ))  

    RHO  = np.zeros((NTX, NTY, NTZ))  
    U = np.zeros((NTX, NTY, NTZ)) 
    V = np.zeros((NTX, NTY, NTZ))  
    W = np.zeros((NTX, NTY, NTZ))  
    P = np.zeros((NTX, NTY, NTZ))  
    TEMP = np.zeros((NTX, NTY, NTZ))  

    
    RHO0 = np.zeros((NTX, NTY, NTZ))  
    U0   = np.zeros((NTX, NTY, NTZ)) 
    V0   = np.zeros((NTX, NTY, NTZ))  
    W0   = np.zeros((NTX, NTY, NTZ))  
    P0   = np.zeros((NTX, NTY, NTZ))  
    TEMP0 = np.zeros((NTX, NTY, NTZ))  


    # mesh (including halo points)

    xf = np.linspace(-order*dxx, Lx+order*dxx, NTX+1)
    yf = np.linspace(-order*dyy, Ly+order*dyy, NTY+1)
    zf = np.linspace(-order*dzz, Lz+order*dzz, NTZ+1)
    
    xf[0] = -order*dxx
    for i in range(1, NTX+1):
        xf[i] = xf[i-1] + dxx
    yf[0] = -order*dyy
    for j in range(1, NTY+1):
        yf[i] = yf[i-1] + dyy
    zf[0] = -order*dzz
    for k in range(1, NTZ+1):
        zf[i] = zf[i-1] + dzz

    xc = 0.5 * (xf[:-1] + xf[1:])
    yc = 0.5 * (yf[:-1] + yf[1:])
    zc = 0.5 * (zf[:-1] + zf[1:])

    XF, YF, ZF = np.meshgrid(xf, yf, zf, indexing='ij')  # FACE 
    XC, YC, ZC = np.meshgrid(xc, yc, zc, indexing='ij')  # MESH

    #print( " Initialise arrays FV or FD") 
    if (FiniteVolume==1):
        for i in range(0, NTX):
            for j in range(0, NTY):
                for k in range(0, NTZ):
                    x1 = XC[i,j,k]
                    y1 = YC[i,j,k]
                    z1 = ZC[i,j,k]
                    RHO[i,j,k]  =  rhoI_func(x1,y1,z1,dxx,dyy,dzz)*o_dvol
                    U[i,j,k]    =  uI_func(x1,y1,z1,dxx,dyy,dzz)*o_dvol
                    V[i,j,k]    =  vI_func(x1,y1,z1,dxx,dyy,dzz)*o_dvol
                    #W[i,j,k]   =  wI_func(x1,y1,z1,dxx,dyy,dzz)*o_dvol
                    P[i,j,k]    =  pI_func(x1,y1,z1,dxx,dyy,dzz)*o_dvol
                    RHO0[i,j,k] = RHO[i,j,k] 
                    U0[i,j,k]   = U[i,j,k] 
                    V0[i,j,k]   = V[i,j,k] 
                    P0[i,j,k]   = P[i,j,k] 
    else:
        RHO  = rho_func(XC,YC,ZC)
        U    = u_func(XC,YC,ZC)
        V    = v_func(XC,YC,ZC)
        #W   = w_func(XC,YC,ZC)
        P    = p_func(XC,YC,ZC)

        RHO0 = rho_func(XC,YC,ZC)
        U0   = u_func(XC,YC,ZC)
        V0   = v_func(XC,YC,ZC)
        #W0   = w_func(XC,YC,ZC)
        P0   = p_func(XC,YC,ZC)

    # Temp
    TEMP  = (P/RHO)*o_Rgas
    TEMP0 = (P0/RHO0)*o_Rgas

    # print( " change to conservative ")
    RHOU = RHO*U
    RHOV = RHO*V
    RHOW = RHO*W
    RHOE = P/(gamma-1) + 0.5*RHO*(U*U + V*V + W*W)

    RHOU0 = RHO*U
    RHOV0 = RHO*V
    RHOW0 = RHO*W
    RHOE0 = P/(gamma-1) + 0.5*RHO*(U*U + V*V + W*W)


    # FLUXES and RHS
    FLX_X   = np.zeros((NTX, NTY, NTZ,NCOMP))  
    FLX_Y   = np.zeros((NTX, NTY, NTZ,NCOMP)) 
    FLX_Z   = np.zeros((NTX, NTY, NTZ,NCOMP))  

    # FLUXES and RHS
    FLX_FAC_X   = np.zeros((NTX+1, NTY, NTZ,NCOMP))  
    FLX_FAC_Y   = np.zeros((NTX, NTY+1, NTZ,NCOMP)) 
    FLX_FAC_Z   = np.zeros((NTX, NTY, NTZ+1,NCOMP))  

    RHS     = np.zeros((NTX, NTY, NTZ,NCOMP))  

    if (EulerFlux ==1):
        FLX_X[:,:,:,0] = RHOU
        FLX_X[:,:,:,1] = RHO*U*U + P
        FLX_X[:,:,:,2] = RHO*U*V
        FLX_X[:,:,:,3] = RHO*U*W
        FLX_X[:,:,:,4] = (RHOE + P)*U

        FLX_Y[:,:,:,0] = RHOV
        FLX_Y[:,:,:,1] = RHO*U*V
        FLX_Y[:,:,:,2] = RHO*V*V + P
        FLX_Y[:,:,:,3] = RHO*V*W
        FLX_Y[:,:,:,4] = (RHOE + P)*V

        FLX_Z[:,:,:,0] = RHOW
        FLX_Z[:,:,:,1] = RHO*U*W
        FLX_Z[:,:,:,2] = RHO*V*W 
        FLX_Z[:,:,:,3] = RHO*W*W + P
        FLX_Z[:,:,:,4] = (RHOE + P)*W

    dt = 1e-6
    CFL  = (u_max+csound)/(dxx/dt)
    CFLv = (visc/rho_min)/(dxx*dxx/dt) 
    print(' CFL=', CFL, ' CFL visc = ',CFLv)
#------------------------------------------------------
    if (Model==0):
        print(" Using central FD   RHS = nabla(F + Fvisc) ")
        if (NavierStokes==1):
            visco_solver.calc_viscfluxes(FiniteVolume,index_start_bc,index_end_bc, visc,cond,o_dx,o_dy,o_dz, U, V, W, TEMP, FLX_X,FLX_Y,FLX_Z)
        # Compute Flux derivatives + Source
        for i in range(ist, ien+1):
            for j in range(jst, jen+1):
                for k in range(kst, ken+1):
                    for nv in range(0,5):
                        dfdx  = derivatives.dfdx_cc(i,j,k,0,FLX_X[:,:,:,nv])*o_dx
                        dfdx += derivatives.dfdx_cc(i,j,k,1,FLX_Y[:,:,:,nv])*o_dy
                        dfdx += derivatives.dfdx_cc(i,j,k,2,FLX_Z[:,:,:,nv])*o_dz
                        RHS[i,j,k,nv] = -dfdx
                    x1 = XC[i,j,k]
                    y1 = YC[i,j,k]
                    z1 = ZC[i,j,k]
                    RHS[i,j,k,0] += S1_func(x1,y1,z1)
                    RHS[i,j,k,1] += S2_func(x1,y1,z1)
                    RHS[i,j,k,2] += S3_func(x1,y1,z1)
                    RHS[i,j,k,3] += S4_func(x1,y1,z1)
                    RHS[i,j,k,4] += S5_func(x1,y1,z1)    
    #------------------------------------------------------
    if (Model==1):
        print(" Using conservative FD   RHS = nabla(F) +nabla(Fvisc)")
        # computer grad(fluxes)
        for i in range(ist, ien+1):
            for j in range(jst, jen+1):
                for k in range(kst, ken+1):
                    for nv in range(0,NCOMP):
                        fl  = derivatives.interp_fd(i  ,j,k,0,nv,FLX_X)
                        fr  = derivatives.interp_fd(i+1,j,k,0,nv,FLX_X)
                        dfdx =  (fr -fl)*o_dx
                        fl  = derivatives.interp_fd(i,j  ,k,1,nv,FLX_Y)
                        fr  = derivatives.interp_fd(i,j+1,k,1,nv,FLX_Y)
                        dfdx +=  (fr -fl)*o_dy
                        fl  = derivatives.interp_fd(i,j,  k,2,nv,FLX_Z)
                        fr  = derivatives.interp_fd(i,j,k+1,2,nv,FLX_Z)
                        dfdx +=  (fr -fl)*o_dz
                        RHS[i,j,k,nv] = -dfdx    
                    x1 = XC[i,j,k]
                    y1 = YC[i,j,k]
                    z1 = ZC[i,j,k]
                    RHS[i,j,k,0] += S1_func(x1,y1,z1)
                    RHS[i,j,k,1] += S2_func(x1,y1,z1)
                    RHS[i,j,k,2] += S3_func(x1,y1,z1)
                    RHS[i,j,k,3] += S4_func(x1,y1,z1)
                    RHS[i,j,k,4] += S5_func(x1,y1,z1)
        if (NavierStokes==1):
            print("  Model 1 NS NRY")
    #------------------------------------------------------
    if (Model==2):
        print(" Using conservative FV")
        # Compute Navier Stokes
        if (NavierStokes==1):
            # Viscous fluxes
            visco_solver.calc_viscfluxes(FiniteVolume,index_start_bc,index_end_bc, visc,cond,o_dx,o_dy,o_dz, U, V, W, TEMP, FLX_FAC_X,FLX_FAC_Y,FLX_FAC_Z)

        # compute EulerFluxes  X
        for i in range(ist, ien+2):
            for j in range(jst, jen+1):
                for k in range(kst, ken+1):
                    for nv in range(0,NCOMP):
                        FLX_FAC_X[i,j,k,nv]  += derivatives.interp_fv(i,j,k,0,nv,FLX_X)
        # compute EulerFluxes  Y
        for i in range(ist, ien+1):
            for j in range(jst, jen+2):
                for k in range(kst, ken+1):
                    for nv in range(0,NCOMP):
                        FLX_FAC_Y[i,j,k,nv]  += derivatives.interp_fv(i,j,k,1,nv,FLX_Y)
        # compute EulerFluxes  Z
        for i in range(ist, ien+1):
            for j in range(jst, jen+1):
                for k in range(kst, ken+2):
                    for nv in range(0,NCOMP):
                        FLX_FAC_Z[i,j,k,nv]  += derivatives.interp_fv(i,j,k,2,nv,FLX_Z)
                        
        # RHS += nabla F)
        for i in range(ist, ien+1):
            for j in range(jst, jen+1):
                for k in range(kst, ken+1):
                    for nv in range(0,NCOMP):
                        RHS[i,j,k,nv] =   (FLX_FAC_X[i,j,k,nv] - FLX_FAC_X[i+1,j,k,nv])*o_dx
                        RHS[i,j,k,nv] +=  (FLX_FAC_Y[i,j,k,nv] - FLX_FAC_Y[i,j+1,k,nv])*o_dy
                        RHS[i,j,k,nv] +=  (FLX_FAC_Z[i,j,k,nv] - FLX_FAC_Z[i,j,k+1,nv])*o_dz
                               
                    RHS[i,j,k,0] += S1I_func(XC[i,j,k],YC[i,j,k],ZC[i,j,k],dxx,dyy,dzz)*o_dvol
                    RHS[i,j,k,1] += S2I_func(XC[i,j,k],YC[i,j,k],ZC[i,j,k],dxx,dyy,dzz)*o_dvol
                    RHS[i,j,k,2] += S3I_func(XC[i,j,k],YC[i,j,k],ZC[i,j,k],dxx,dyy,dzz)*o_dvol
                    RHS[i,j,k,3] += S4I_func(XC[i,j,k],YC[i,j,k],ZC[i,j,k],dxx,dyy,dzz)*o_dvol
                    RHS[i,j,k,4] += S5I_func(XC[i,j,k],YC[i,j,k],ZC[i,j,k],dxx,dyy,dzz)*o_dvol

        
    #------------------------------------------------------
    # Advance Solution (All methods the same)
    RHO  += dt*RHS[:,:,:,0] 
    RHOU += dt*RHS[:,:,:,1] 
    RHOV += dt*RHS[:,:,:,2] 
    RHOW += dt*RHS[:,:,:,3] 
    RHOE += dt*RHS[:,:,:,4] 
        
    # from U-->Q
    for i in range(ist, ien+1):
        for j in range(jst, jen+1):
            for k in range(kst, ken+1):
                U[i,j,k] = RHOU[i,j,k]/RHO[i,j,k]
                V[i,j,k] = RHOV[i,j,k]/RHO[i,j,k]
                W[i,j,k] = RHOW[i,j,k]/RHO[i,j,k]
                kin = U[i,j,k]*U[i,j,k]+V[i,j,k]*V[i,j,k] + W[i,j,k]*W[i,j,k]
                rhoe_int = RHOE[i,j,k] - RHO[i,j,k]*kin*0.5
                P[i,j,k] = rhoe_int*(gamma-1)

    # cal errors
    NORME = 1e6
    if (ErrorQ==1):
        ERR_RHO   = RHO - RHO0
        ERR_U     = U   - U0
        ERR_V     = V   - V0
        ERR_W     = W   - W0
        ERR_P     = P   - P0
    else:
        ERR_RHO   = RHO - RHO0
        ERR_U     = RHOU   - RHOU0
        ERR_V     = RHOV   - RHOV0
        ERR_W     = RHOW   - RHOW0
        ERR_P     = RHOE   - RHOE0

    
    # error only in internal points
    HALO2 = 2*HALO
    interior1 = ERR_RHO[HALO2:-HALO2, HALO2:-HALO2, HALO2:-HALO2] 
    interior2 = ERR_U[HALO2:-HALO2, HALO2:-HALO2, HALO2:-HALO2] 
    interior3 = ERR_V[HALO2:-HALO2, HALO2:-HALO2, HALO2:-HALO2] 
    interior4 = ERR_W[HALO2:-HALO2, HALO2:-HALO2, HALO2:-HALO2] 
    interior5 = ERR_P[HALO2:-HALO2, HALO2:-HALO2, HALO2:-HALO2] 
    
    if  (ErrorMeasure==1):
        lerror_rho = np.mean( np.abs(interior1) )
        lerror_u   = np.mean( np.abs(interior2)) 
        lerror_v   = np.mean( np.abs(interior3) ) 
        lerror_w   = np.mean( np.abs(interior4) ) 
        lerror_p   = np.mean( np.abs(interior5))  

    elif (ErrorMeasure==0):
        lerror_rho = np.max( np.abs(interior1)    )
        lerror_u   = np.max( np.abs(interior2)   ) 
        lerror_v   = np.max( np.abs(interior3)   ) 
        lerror_w   = np.max( np.abs(interior4)   ) 
        lerror_p   = np.max( np.abs(interior5)   ) 
    
    elif (ErrorMeasure==2):
        lerror_rho = np.mean(interior1**2 )**0.5
        lerror_u   = np.mean(interior2**2 )**0.5
        lerror_v   = np.mean(interior3**2 )**0.5
        lerror_w   = np.mean(interior4**2 )**0.5
        lerror_p   = np.mean(interior5**2 )**0.5
    else:
        print("not error measure specified")

    
    errors_rho.append(lerror_rho)
    errors_u.append(lerror_u)
    errors_v.append(lerror_v)
    errors_w.append(lerror_w)
    errors_p.append(lerror_p)

    print(" N =",N, "---- ErrorMeasure   = ", ErrorMeasure)
   #print(" ERROR rho= ",lerror_rho)
    print(" ERROR u= ",lerror_u)
    # print(" ERROR v= ",lerror_v)
    # print(" ERROR w= ",lerror_w)
    # print(" ERROR p= ",lerror_p)


# Plot convergence
print(" ... PLOTTING .... ")

plt.figure(figsize=(8, 6))
e0    = errors_rho[0]
e0u   = errors_u[0]
e0v   = errors_v[0]
e0w   = errors_w[0]
e0p   = errors_p[0]


plt.loglog(resolutions, errors_rho/e0, '-o', color='green',label='rho',markersize=10)
plt.loglog(resolutions,  errors_u/e0u, '>', color='red',label='  u',markersize=10)
#plt.loglog(resolutions,  errors_v/e0v, '*', color='cyan',label='  v')
#plt.loglog(resolutions,  errors_w/e0w, '+', color='red',label='  w')
plt.loglog(resolutions,  errors_p/e0p, 'D', color='blue',label='  p',markersize=10)


# Add ideal convergence lines
resolutions = np.array(resolutions)  # Convert list to NumPy array
N_ref = resolutions[0]
ideal_1st = (resolutions / N_ref) ** -1
ideal_2nd = (resolutions / N_ref) ** -2
ideal_4th = (resolutions / N_ref) ** -4
ideal_6th = (resolutions / N_ref) ** -6

plt.loglog(resolutions, ideal_1st, linestyle='--', color='gray')
plt.loglog(resolutions, ideal_2nd, linestyle='--', color='gray')
plt.loglog(resolutions, ideal_4th, linestyle='--', color='gray')
plt.loglog(resolutions, ideal_6th, linestyle='--', color='gray')


# Add text labels at the last point of each line
x_text = resolutions[-1]*0.8
plt.text(x_text, ideal_1st[-1], "1st", fontsize=10, va='bottom', ha='left')
plt.text(x_text, ideal_2nd[-1], "2nd", fontsize=10, va='bottom', ha='left')
plt.text(x_text, ideal_4th[-1], "4th", fontsize=10, va='bottom', ha='left')
plt.text(x_text, ideal_6th[-1], "6th", fontsize=10, va='bottom', ha='left')



ax = plt.gca()  # get current axis
ax.xaxis.set_major_formatter(ticker.ScalarFormatter())
ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%d'))
#ax.ticklabel_format(style='plain', axis='x')  # ensure no scientific notation

#plt.grid(True, which='both', ls='--')
plt.xlabel('Resolution (N)')
plt.ylabel('Normalised  Error')
#plt.title('Error Convergence')
plt.legend()

# Optional: fit and plot convergence rate
p = np.polyfit(np.log(resolutions), np.log(errors_rho), 1)
order = -p[0]      

q = np.polyfit(np.log(resolutions), np.log(errors_u), 1)
order_u = -q[0]

r = np.polyfit(np.log(resolutions), np.log(errors_p), 1)
order_p = -r[0]


plt.text(resolutions[0], 0.01*errors_u[0]/e0u, f"Order ≈ {order_u:.2f}", fontsize=12)
plt.show()



print(" order(rho)= ",order, " order(u)= ",order_u," order(p)= ",order_p)

# k = 10  # fixed z-index
# plt.imshow(U[:, :, k], origin='lower', cmap='viridis')
# plt.colorbar(label='U')
# plt.title(f"Slice at k = {k}")
# plt.xlabel('j')
# plt.ylabel('i')
# plt.show()
