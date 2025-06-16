import sympy as sp
import numpy as np
import matplotlib.pyplot as plt
import derivatives


NavierStokes = 0  # if 1 solve Navier Stokes
Model        = 0  # 0:central  1:conservative flux/FD  2:conservative flux/FV
FiniteVolume = 0  # 0:FD      1: FV
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
w_oo   = 0.0

csound = u_oo/Ma
p_oo   = rho_oo*csound*csound/gamma
T_oo   = p_oo/(Rgas*rho_oo)

visc = rho_oo*u_oo*L_oo/Re

cond = visc*Cp/Pr

# print(" csound = ", csound)
# print(" visc = ", visc, " cond= ", cond)
# print(" p_oo = ", p_oo, " T_oo = ",T_oo)
# print(" Mach      = ", u_oo/csound)
# print(" Reynolds  = ", rho_oo*u_oo*L_oo/visc)



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


# ds = 0.01
# x1 = 0.5
# y1 = 0.5
# z1 = 0.5
# dvol = ds*ds*ds

# rho_val = rhoI.evalf(subs={x: x1, y: y1, z: z1, dx: ds, dy: ds, dz: ds})
# print("rho FINITE VOLUME     rhof  = 1/V int rho dV:", rho_val/dvol)
# print("rho FINITE DIFFERENCE rho   = rho(x,y,z    ):", rho_func(x1,y1,z1))
# print("rho FINITE VOLUME     rhoI_func  :", rhoI_func(x1,y1,z1,ds,ds,ds)/dvol)


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


if NavierStokes  >0:
    # Visc Fluxes x
    F1viscx = 0
    F2viscx = tau_xx
    F3viscx = tau_xy
    F4viscx = tau_xz
    F5viscx = cond*sp.diff(T, x) + u*tau_xx + v*tau_xy + w*tau_xz
    # Visc Fluxes y
    F1viscy = 0
    F2viscy = tau_xy
    F3viscy = tau_yy
    F4viscy = tau_yz
    F5viscy = cond*sp.diff(T, y) + u*tau_xy + v*tau_yy + w*tau_yz
    # Visc Fluxes z
    F1viscz = 0
    F2viscz = tau_xz
    F3viscz = tau_yz
    F4viscz = tau_zz
    F5viscz = cond*sp.diff(T, z) + u*tau_xz + v*tau_yz + w*tau_zz
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

# functions
S1_func   = sp.lambdify((x, y, z), S1, modules=["numpy"])
S2_func   = sp.lambdify((x, y, z), S2, modules=["numpy"])
S3_func   = sp.lambdify((x, y, z), S3, modules=["numpy"])
S4_func   = sp.lambdify((x, y, z), S4, modules=["numpy"])
S5_func   = sp.lambdify((x, y, z), S5, modules=["numpy"])

# Compute the integral
if (FiniteVolume==1):
    S1I   = sp.integrate(S1  , z_bounds, y_bounds, x_bounds)
    S2I   = sp.integrate(S2  , z_bounds, y_bounds, x_bounds)
    S3I   = sp.integrate(S3  , z_bounds, y_bounds, x_bounds)
    S4I   = sp.integrate(S4  , z_bounds, y_bounds, x_bounds)
    S5I   = sp.integrate(S5  , z_bounds, y_bounds, x_bounds)

    S1I_func = sp.lambdify((x, y, z, dx, dy, dz), S1I, modules=["numpy"])
    S2I_func = sp.lambdify((x, y, z, dx, dy, dz), S2I, modules=["numpy"])
    S3I_func = sp.lambdify((x, y, z, dx, dy, dz), S3I, modules=["numpy"])
    S4I_func = sp.lambdify((x, y, z, dx, dy, dz), S4I, modules=["numpy"])
    S5I_func = sp.lambdify((x, y, z, dx, dy, dz), S5I, modules=["numpy"])

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
#resolutions = [8, 16, 32, 64]
resolutions = [8, 16, 32]
#resolutions = [4]

errors_rho  = []
errors_u    = []
errors_v    = []
errors_w    = []
errors_p    = []

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
    ist = HALO          # 1/2/3 for 2/4/6 (first inner point)
    ien = ist + Nx - 1  # 8/9/10
    jst = HALO         
    jen = ist + Ny - 1
    kst = HALO          
    ken = kst + Nz - 1

    NTX = Nx + 2*HALO 
    NTY = Ny + 2*HALO
    NTZ = Nz + 2*HALO

    # print( "ist jst kst ", ist, jst, kst)
    # print( "ien jen ken ", ien, jen, ken)
    # print( "NTX NTY NTZ ", NTX, NTY, NTZ)

    # 3D array with all elements initialized to 0.0
    RHOU = np.zeros((NTX, NTY, NTZ))  
    RHOV = np.zeros((NTX, NTY, NTZ))  
    RHOW = np.zeros((NTX, NTY, NTZ))  
    RHOET = np.zeros((NTX, NTY, NTZ))  

    RHS_RHO  = np.zeros((NTX, NTY, NTZ))  
    RHS_RHOU = np.zeros((NTX, NTY, NTZ))  
    RHS_RHOV = np.zeros((NTX, NTY, NTZ))  
    RHS_RHOW = np.zeros((NTX, NTY, NTZ))  
    RHS_RHOET = np.zeros((NTX, NTY, NTZ))  

    RHO  = np.zeros((NTX, NTY, NTZ))  
    U = np.zeros((NTX, NTY, NTZ)) 
    V = np.zeros((NTX, NTY, NTZ))  
    W = np.zeros((NTX, NTY, NTZ))  
    P = np.zeros((NTX, NTY, NTZ))  
    
    RHO0 = np.zeros((NTX, NTY, NTZ))  
    U0   = np.zeros((NTX, NTY, NTZ)) 
    V0   = np.zeros((NTX, NTY, NTZ))  
    W0   = np.zeros((NTX, NTY, NTZ))  
    P0   = np.zeros((NTX, NTY, NTZ))  

    # mesh (including halo points)

    xf = np.linspace(-HALO*dxx, Lx+HALO*dxx, NTX+1)
    yf = np.linspace(-HALO*dyy, Ly+HALO*dyy, NTY+1)
    zf = np.linspace(-HALO*dzz, Lz+HALO*dzz, NTZ+1)
    
    xf[0] = -HALO*dxx
    for i in range(1, NTX+1):
        xf[i] = xf[i-1] + dxx
    yf[0] = -HALO*dyy
    for j in range(1, NTY+1):
        yf[i] = yf[i-1] + dyy
    zf[0] = -HALO*dzz
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
        #W   = w_func(XC,YC,ZC)
        P0   = p_func(XC,YC,ZC)


    # print( " change to conservative ")
    RHOU = RHO*U
    RHOV = RHO*V
    RHOW = RHO*W
    RHOE = P/(gamma-1) + 0.5*RHO*(U*U + V*V + W*W)

    # auxiliar variables

    RHOUU = np.zeros((NTX, NTY, NTZ))  
    RHOUU = RHO*U*U + P
    RHOUV = np.zeros((NTX, NTY, NTZ))  
    RHOUV = RHO*U*V
    RHOUW = np.zeros((NTX, NTY, NTZ))  
    RHOUW = RHO*U*W
    RHOVV = np.zeros((NTX, NTY, NTZ))  
    RHOVV = RHO*V*V + P
    RHOVW = np.zeros((NTX, NTY, NTZ))  
    RHOVW = RHO*V*W
    RHOWW = np.zeros((NTX, NTY, NTZ))  
    RHOWW = RHO*W*W + P

    FEX = np.zeros((NTX, NTY, NTZ))  
    FEX = (RHOE + P)*U
    FEY = np.zeros((NTX, NTY, NTZ))  
    FEY = (RHOE + P)*V
    FEZ = np.zeros((NTX, NTY, NTZ))  
    FEZ = (RHOE + P)*W 

    # FLUXES
    FLX_X   = np.zeros((NTX, NTY, NTZ,5))  
    FLX_Y   = np.zeros((NTX, NTY, NTZ,5)) 
    FLX_Z   = np.zeros((NTX, NTY, NTZ,5))  
    RHS     = np.zeros((NTX, NTY, NTZ,5))  

    # 
    for i in range(0, NTX):
        for j in range(0, NTY):
            for k in range(0, NTZ):
                FLX_X[i,j,k,0] = RHOU[i,j,k]
                FLX_X[i,j,k,1] = RHOUU[i,j,k]
                FLX_X[i,j,k,2] = RHOUV[i,j,k]
                FLX_X[i,j,k,3] = RHOUW[i,j,k]
                FLX_X[i,j,k,4] = FEX[i,j,k]
                FLX_Y[i,j,k,0] = RHOV[i,j,k]
                FLX_Y[i,j,k,1] = RHOUV[i,j,k]
                FLX_Y[i,j,k,2] = RHOVV[i,j,k]
                FLX_Y[i,j,k,3] = RHOVW[i,j,k]
                FLX_Y[i,j,k,4] = FEY[i,j,k]
                FLX_Z[i,j,k,0] = RHOW[i,j,k]
                FLX_Z[i,j,k,1] = RHOUW[i,j,k]
                FLX_Z[i,j,k,2] = RHOVW[i,j,k]
                FLX_Z[i,j,k,3] = RHOWW[i,j,k]
                FLX_Z[i,j,k,4] = FEZ[i,j,k]

    # CENTRAL DERIVATIVE

    dt = 0.00001
    CFL = (u_max+csound)/(dxx/dt)
    print(' CFL=', CFL)
#------------------------------------------------------
    if (Model==0):
        print(" Using central FD")
        # Compute Euler
        for i in range(ist, ien+1):
            for j in range(jst, jen+1):
                for k in range(kst, ken+1):
                    # rho
                    dfdx  = derivatives.dfdx_cc(i,j,k,0,RHOU)*o_dx
                    dfdx += derivatives.dfdx_cc(i,j,k,1,RHOV)*o_dy
                    dfdx += derivatives.dfdx_cc(i,j,k,2,RHOW)*o_dz
                    source = S1_func(XC[i,j,k],YC[i,j,k],ZC[i,j,k])
                    RHS_RHO[i,j,k] = -dfdx + source
                    RHO[i,j,k]  += dt*RHS_RHO[i,j,k]
                    # rhou
                    dfdx  = derivatives.dfdx_cc(i,j,k,0,RHOUU)*o_dx
                    dfdx += derivatives.dfdx_cc(i,j,k,1,RHOUV)*o_dy
                    dfdx += derivatives.dfdx_cc(i,j,k,2,RHOUW)*o_dz
                    source = S2_func(XC[i,j,k],YC[i,j,k],ZC[i,j,k])
                    RHS_RHOU[i,j,k] = -dfdx + source
                    RHOU[i,j,k]  += dt*RHS_RHOU[i,j,k]
                    # rhov
                    dfdx  = derivatives.dfdx_cc(i,j,k,0,RHOUV)*o_dx
                    dfdx += derivatives.dfdx_cc(i,j,k,1,RHOVV)*o_dy
                    dfdx += derivatives.dfdx_cc(i,j,k,2,RHOVW)*o_dz
                    source = S3_func(XC[i,j,k],YC[i,j,k],ZC[i,j,k])
                    RHS_RHOV[i,j,k] = -dfdx + source
                    RHOV[i,j,k]  += dt*RHS_RHOV[i,j,k]
                    # rhow
                    dfdx  = derivatives.dfdx_cc(i,j,k,0,RHOUW)*o_dx
                    dfdx += derivatives.dfdx_cc(i,j,k,1,RHOVW)*o_dy
                    dfdx += derivatives.dfdx_cc(i,j,k,2,RHOWW)*o_dz
                    source = S4_func(XC[i,j,k],YC[i,j,k],ZC[i,j,k])
                    RHS_RHOW[i,j,k] = -dfdx + source
                    RHOW[i,j,k]  += dt*RHS_RHOW[i,j,k]
                    
                    # rhoe
                    dfdx  = derivatives.dfdx_cc(i,j,k,0,FEX)*o_dx
                    dfdx += derivatives.dfdx_cc(i,j,k,1,FEY)*o_dy
                    dfdx += derivatives.dfdx_cc(i,j,k,2,FEZ)*o_dz
                    source = S5_func(XC[i,j,k],YC[i,j,k],ZC[i,j,k])
                    RHS_RHOET[i,j,k] = -dfdx + source
                    RHOE[i,j,k]  += dt*RHS_RHOET[i,j,k]
        # Compute Navier Stokes
        if (NavierStokes==1):
            # Viscous
            print(" solve N-S")
        
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

        # errors
        # RHS_RHO[i,j,k]  = RHO[i,j,k] - RHO0[i,j,k]
        # RHS_RHOU[i,j,k] = U[i,j,k]   - U0[i,j,k]
        # RHS_RHOV[i,j,k] = V[i,j,k]   - V0[i,j,k]
        # RHS_RHOW[i,j,k] = W[i,j,k] - W0[i,j,k]
        # RHS_RHOET[i,j,k] = P[i,j,k]   - P0[i,j,k]

        RHS_RHO   = RHO - RHO0
        RHS_RHOU  = U   - U0
        RHS_RHOV  = V   - V0
        RHS_RHOW  = W   - W0
        RHS_RHOET = P   - P0
#------------------------------------------------------
    if (Model==1):
        print(" Using conservative FD")
        for i in range(ist, ien+1):
            for j in range(jst, jen+1):
                for k in range(kst, ken+1):
                    # loop over eqns
                    for nv in range(0,5):
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
                  
                    RHS_RHO[i,j,k] = RHS[i,j,k,0]  + S1_func(XC[i,j,k],YC[i,j,k],ZC[i,j,k])
                    RHO[i,j,k] +=  dt * RHS_RHO[i,j,k]

                    RHS_RHOU[i,j,k] = RHS[i,j,k,1] + S2_func(XC[i,j,k],YC[i,j,k],ZC[i,j,k])
                    RHOU[i,j,k]  += dt*RHS_RHOU[i,j,k]

                    RHS_RHOV[i,j,k] = RHS[i,j,k,2] + S3_func(XC[i,j,k],YC[i,j,k],ZC[i,j,k])
                    RHOV[i,j,k]  += dt*RHS_RHOV[i,j,k]

                    RHS_RHOW[i,j,k] = RHS[i,j,k,3] + S4_func(XC[i,j,k],YC[i,j,k],ZC[i,j,k])
                    RHOW[i,j,k]  += dt*RHS_RHOW[i,j,k]

                    RHS_RHOET[i,j,k] = RHS[i,j,k,4]+ S5_func(XC[i,j,k],YC[i,j,k],ZC[i,j,k])
                    RHOE[i,j,k]  += dt*RHS_RHOET[i,j,k]

                    # errores
                    RHS_RHO[i,j,k]  = RHO[i,j,k] - RHO0[i,j,k]
                    U[i,j,k] = RHOU[i,j,k]/RHO[i,j,k]
                    RHS_RHOU[i,j,k] = U[i,j,k]   - U0[i,j,k]
                    V[i,j,k] = RHOV[i,j,k]/RHO[i,j,k]
                    RHS_RHOV[i,j,k] = V[i,j,k]   - V0[i,j,k]
                    W[i,j,k] = RHOW[i,j,k]/RHO[i,j,k]
                    RHS_RHOW[i,j,k] = W[i,j,k]   - W0[i,j,k]
                    kin = U[i,j,k]*U[i,j,k]+V[i,j,k]*V[i,j,k] + W[i,j,k]*W[i,j,k]
                    rhoe_int = RHOE[i,j,k] - RHO[i,j,k]*kin*0.5
                    P[i,j,k] = rhoe_int*(gamma-1)
                    RHS_RHOET[i,j,k] = P[i,j,k]   - P0[i,j,k]
    #------------------------------------------------------
    if (Model==2):
        print(" Using conservative FV")
        for i in range(ist, ien+1):
            for j in range(jst, jen+1):
                for k in range(kst, ken+1):
                    # loop over eqns
                    for nv in range(0,5):
                        fl  = derivatives.interp_fv(i  ,j,k,0,nv,FLX_X)
                        fr  = derivatives.interp_fv(i+1,j,k,0,nv,FLX_X)
                        dfdx =  (fr -fl)*o_dx
                        fl  = derivatives.interp_fv(i,j  ,k,1,nv,FLX_Y)
                        fr  = derivatives.interp_fv(i,j+1,k,1,nv,FLX_Y)
                        dfdx +=  (fr -fl)*o_dy
                        fl  = derivatives.interp_fv(i,j,  k,2,nv,FLX_Z)
                        fr  = derivatives.interp_fv(i,j,k+1,2,nv,FLX_Z)
                        dfdx +=  (fr -fl)*o_dz
                        RHS[i,j,k,nv] = -dfdx
                  
                    RHS_RHO[i,j,k] = RHS[i,j,k,0]  + S1I_func(XC[i,j,k],YC[i,j,k],ZC[i,j,k],dxx,dyy,dzz)*o_dvol
                    RHO[i,j,k] +=  dt * RHS_RHO[i,j,k]

                    RHS_RHOU[i,j,k] = RHS[i,j,k,1] + S2I_func(XC[i,j,k],YC[i,j,k],ZC[i,j,k],dxx,dyy,dzz)*o_dvol
                    RHOU[i,j,k]  += dt*RHS_RHOU[i,j,k]

                    RHS_RHOV[i,j,k] = RHS[i,j,k,2] + S3I_func(XC[i,j,k],YC[i,j,k],ZC[i,j,k],dxx,dyy,dzz)*o_dvol
                    RHOV[i,j,k]  += dt*RHS_RHOV[i,j,k]

                    RHS_RHOW[i,j,k] = RHS[i,j,k,3] + S4I_func(XC[i,j,k],YC[i,j,k],ZC[i,j,k],dxx,dyy,dzz)*o_dvol
                    RHOW[i,j,k]  += dt*RHS_RHOW[i,j,k]

                    RHS_RHOET[i,j,k] = RHS[i,j,k,4]+ S5I_func(XC[i,j,k],YC[i,j,k],ZC[i,j,k],dxx,dyy,dzz)*o_dvol
                    RHOE[i,j,k]  += dt*RHS_RHOET[i,j,k]

                    # errores
                    RHS_RHO[i,j,k]  = RHO[i,j,k] - RHO0[i,j,k]
                    U[i,j,k] = RHOU[i,j,k]/RHO[i,j,k]
                    RHS_RHOU[i,j,k] = U[i,j,k]   - U0[i,j,k]
                    V[i,j,k] = RHOV[i,j,k]/RHO[i,j,k]
                    RHS_RHOV[i,j,k] = V[i,j,k]   - V0[i,j,k]
                    W[i,j,k] = RHOW[i,j,k]/RHO[i,j,k]
                    RHS_RHOW[i,j,k] = W[i,j,k]   - W0[i,j,k]
                    kin = U[i,j,k]*U[i,j,k]+V[i,j,k]*V[i,j,k] + W[i,j,k]*W[i,j,k]
                    rhoe_int = RHOE[i,j,k] - RHO[i,j,k]*kin*0.5
                    P[i,j,k] = rhoe_int*(gamma-1)
                    RHS_RHOET[i,j,k] = P[i,j,k]   - P0[i,j,k]

    #------------------------------------------------------

    l2error_rho = np.mean( RHS_RHO*RHS_RHO     )
    l2error_u   = np.mean( RHS_RHOU*RHS_RHOU   ) 
    l2error_v   = np.mean( RHS_RHOV*RHS_RHOV   ) 
    l2error_w   = np.mean( RHS_RHOW*RHS_RHOW   ) 
    l2error_p   = np.mean( RHS_RHOET*RHS_RHOET ) 

    #  max errors
    l2error_rho = np.max( np.abs(RHS_RHO)     )
    l2error_u   = np.max( np.abs(RHS_RHOU)    ) 
    l2error_p   = np.max( np.abs(RHS_RHOET)   ) 

    
    errors_rho.append(l2error_rho)
    errors_u.append(l2error_u)
    errors_v.append(l2error_v)
    errors_w.append(l2error_w)
    errors_p.append(l2error_p)

    print(" N =",N, "----")
    print(" L2 ERROR rho= ",l2error_rho)
    print(" L2 ERROR u= ",l2error_u)
    print(" L2 ERROR v= ",l2error_v)
    print(" L2 ERROR w= ",l2error_w)
    print(" L2 ERROR p= ",l2error_p)


# Plot convergence
print(" ... PLOTTING .... ")

plt.figure(figsize=(8, 6))
e0    = errors_rho[0]
e0u   = errors_u[0]
e0v   = errors_v[0]
e0w   = errors_w[0]
e0p   = errors_p[0]



plt.loglog(resolutions, errors_rho/e0, '-o', color='green',label='rho')
plt.loglog(resolutions,  errors_u/e0u, '>', color='red',label='  u')
#plt.loglog(resolutions,  errors_v/e0v, '*', color='cyan',label='  v')
#plt.loglog(resolutions,  errors_w/e0w, '+', color='red',label='  w')
plt.loglog(resolutions,  errors_p/e0p, 'D', color='blue',label='  p')


# Add ideal convergence lines
resolutions = np.array(resolutions)  # Convert list to NumPy array
N_ref = resolutions[0]
ideal_2nd = (resolutions / N_ref) ** -2
ideal_4th = (resolutions / N_ref) ** -4
ideal_6th = (resolutions / N_ref) ** -6

plt.loglog(resolutions, ideal_2nd, linestyle='--', color='gray', label='2nd-order')
plt.loglog(resolutions, ideal_4th, linestyle='--', color='gray', label='4th-order')
plt.loglog(resolutions, ideal_6th, linestyle='--', color='gray', label='6th-order')


plt.grid(True, which='both', ls='--')
plt.xlabel('Resolution (N)')
plt.ylabel('Normalised L2 Error')
plt.title('Error Convergence')
plt.legend()

# Optional: fit and plot convergence rate
p = np.polyfit(np.log(resolutions), np.log(errors_rho), 1)
order = -p[0]      
plt.text(resolutions[0], errors_rho[0], f"Order ≈ {order:.2f}", fontsize=12)

print(" Order = ",order)

plt.show()