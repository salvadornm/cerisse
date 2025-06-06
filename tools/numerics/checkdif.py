import sympy as sp
import numpy as np
import matplotlib.pyplot as plt

# Symbolic function for exact derivative
x = sp.symbols('x')
f_sym = sp.sin(2 * sp.pi * x)
df_sym = sp.diff(f_sym, x)
f_func_raw = sp.lambdify((x), f_sym, modules=["numpy"])
f_func = lambda val: np.float64(f_func_raw(val))

df_func_raw = sp.lambdify((x), df_sym, modules=["numpy"])
df_func = lambda val: np.float64(df_func_raw(val))

# Symbolic integral (indefinite)
F_sym = sp.integrate(f_sym, x)

# Create a numerical version of the integral
F_func_raw = sp.lambdify((x), F_sym, modules=["numpy"])
F_func = lambda val: np.float64(F_func_raw(val))

# Resolutions to test (finer will not work)
resolutions = [16, 32, 64, 128, 256, 512]
#resolutions = [20]

errors       = []
errors_fc    = []
errors_fv    = []
errors_int   = []
errors_intfv = []


#  derivatives coefficients CELL centred  df(i)   ------------------------

# 2nd order
coef_cc2 = np.array([-0.5, 0.5],dtype=np.float64) 

# 4th order
coef_cc4 = np.array([1.0 / 12 , -2.0 / 3 , 2.0 / 3, -1.0 / 12  ],dtype=np.float64) 

# 6th order
coef_cc6 = np.array([-1.0 / 60 , 3.0/20,  -3.0 / 4 , 3.0 / 4, -3.0/20,  1.0 / 60  ],dtype=np.float64) 

#  derivatives coefficients FACE centred  df(i+1/2) ------------------------

# 2nd order  FD-based and FV-based
coef_fd2 =  np.array([-1, 1],dtype=np.float64) 
coef_fv2 =  np.array([-1, 1],dtype=np.float64) 

# 4th order
coef_fd4 =  np.array([1.0/24, -9.0/8, 9.0/8, -1.0/24],dtype=np.float64) 
coef_fv4 =  np.array([1.0/12, -5.0/4, 5.0/4, -1.0/12],dtype=np.float64) 

# 6th order
coef_fd6 = np.array([-3.0/640, 25.0/384, -75.0/64, 75.0/64, -25.0/384, 3.0/640],dtype=np.float64) 
coef_fv6 = np.array([-1.0/90, 5.0/36, -49.0/36, 49.0/36, -5.0/36, 1.0/90],dtype=np.float64) 

#  interpolation  FD-based
int_fd2 =  np.array([0.5, 0.5], dtype=np.float64) 
int_fd4 =  np.array([-1.0/16, 9.0/16, 9.0/16, -1.0/16], dtype=np.float64) 
int_fd6 =  np.array([3.0/256, -25.0/256, 75.0/128, 75.0/128,-25.0/256, 3.0/256], dtype=np.float64) 

#  interpolation  FV-based
int_fv2 =  np.array([0.5, 0.5], dtype=np.float64) 
int_fv4 =  np.array([-1.0/12, 7.0/12, 7.0/12, -1.0/12], dtype=np.float64) 
int_fv6 =  np.array([1.0/60, -8.0/60, 37.0/60, 37.0/60,-8.0/60, 1.0/60], dtype=np.float64) 

# domain size
L = 1.0

order = 6  #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
HALO = int(order/2)
stencil = 2*HALO
print(" stencil=",stencil, " halo ",HALO)
coef_cc = np.zeros(stencil)   # dfdx coefs at Cell 
coef_fd = np.zeros(stencil)   # dfdx coefs at face based on fd  
coef_fv = np.zeros(stencil)   # dfdx coefs at face based on fv
int_fd = np.zeros(stencil)    # f int coefs at face based on fd
int_fv = np.zeros(stencil)    # f int coefs at face based on fv

if (order==2):
    coef_cc = coef_cc2
    coef_fd = coef_fd2
    coef_fv = coef_fv2
    int_fd = int_fd2
    int_fv = int_fv2
elif (order==4):
    coef_cc = coef_cc4
    coef_fd = coef_fd4
    coef_fv = coef_fv4
    int_fd = int_fd4
    int_fv = int_fv4
elif (order==6):
    coef_cc = coef_cc6
    coef_fd = coef_fd6
    coef_fv = coef_fv6
    int_fd = int_fd6
    int_fv = int_fv6
else:
   raise ValueError("Unsupported order. Only orders 2, 4, and 6 are allowed.")
    
for N in resolutions:
    dx = L / N

    # cells are centred internal point   x[0]= dx/2  there are N cells
    x_cell = np.zeros(N)
    for i in range(0, N):  
        x_cell[i] =  i*dx + 0.5*dx

    # faces ar elocated at   x[0]= 0  x[1] = dx ...  there are N+1 faces
    x_face = np.zeros(N+1)
    for i in range(0, N+1):  
        x_face[i] =  i*dx

    # exact derivative  at cell location
    df_exact = df_func(x_cell)

    # exact derivative in at faces location
    df_exact_face = df_func(x_face)
    
    # function evaluated at cell centred  (Finite Difference)
    f_vals= f_func(x_cell)

    # create FV values
    f_cell = np.zeros(N)
    for i in range(0, N):  
        f_cell[i] =   ( F_func(x_face[i+1]) - F_func(x_face[i]) ) /dx
        

    # function evaluated at faces
    f_exact_face = f_func(x_face)

    # derivative arrays
    df_pred         = np.zeros(N)
    df_pred_face    = np.zeros(N+1)
    df_pred_faceFV  = np.zeros(N+1)

    # interpolation arrays
    f_pred_face     =  np.zeros(N+1)
    f_pred_faceFV   =  np.zeros(N+1)

        
    istart = HALO        # first inner cell point   1/2/3
    iend   = N-HALO -1   # last  inner cell point   N-2/N-3/N-4

    # Interpolation at face-centred   -----------------------
    for i in range(istart, iend + 2):
        for j in range(0, stencil):  # 0,1
            f_pred_face[i]   = f_pred_face[i]   + int_fd[j]*f_vals[i-HALO+j]   #  FD - based
            f_pred_faceFV[i] = f_pred_faceFV[i] + int_fv[j]*f_cell[i-HALO+j]   #  FV - based
    
    # Central derivatives at face-centred   -----------------------
    # derivative at df_face[i] = f[i] -f[i-1] using 2nd order

    for i in range(istart, iend + 2):
        for j in range(0, stencil):  # 0,1
            df_pred_face[i]   = df_pred_face[i]   + coef_fd[j]*f_vals[i-HALO+j]   #  FD - based
            df_pred_faceFV[i] = df_pred_faceFV[i] + coef_fv[j]*f_cell[i-HALO+j]   #  FV - based
   
    # BC (manually could be done better) is periodic
    if (order==2):
        df_pred_face[0]  = coef_fd[1]*f_vals[0] + coef_fd[0]* f_vals[N-1]
        df_pred_face[N]  = df_pred_face[0]

        f_pred_face[0]   = int_fd[1]*f_vals[0] + int_fd[0]* f_vals[N-1]
        f_pred_face[N]   = f_pred_face[0]

        df_pred_faceFV[0]   = coef_fv[1]*f_cell[0] + coef_fv[0]* f_cell[N-1]
        df_pred_faceFV[N]   = df_pred_faceFV[0]

        f_pred_faceFV[0]   = int_fv[1]*f_cell[0] + int_fv[0]* f_cell[N-1]
        f_pred_faceFV[N]   = f_pred_faceFV[0]

    elif (order==4):  
        df_pred_face[0]   = coef_fd[3]*f_vals[1] + coef_fd[2]* f_vals[0] + coef_fd[1]*f_vals[N-1] + coef_fd[0]* f_vals[N-2] 
        df_pred_face[1]   = coef_fd[3]*f_vals[2] + coef_fd[2]* f_vals[1] + coef_fd[1]*f_vals[0]   + coef_fd[0]* f_vals[N-1]
        df_pred_face[N-1] = df_pred_face[1]  
        df_pred_face[N]   = df_pred_face[0]
 
        f_pred_face[0]   = int_fd[3]*f_vals[1] + int_fd[2]* f_vals[0] + int_fd[1]*f_vals[N-1] + int_fd[0]* f_vals[N-2] 
        f_pred_face[1]   = int_fd[3]*f_vals[2] + int_fd[2]* f_vals[1] + int_fd[1]*f_vals[0]   + int_fd[0]* f_vals[N-1]
        f_pred_face[N-1] = int_fd[3]*f_vals[0] + int_fd[2]* f_vals[N-1] + int_fd[1]*f_vals[N-2] + int_fd[0]* f_vals[N-3]  
        f_pred_face[N]   = f_pred_face[0]

        df_pred_faceFV[0]   = coef_fv[3]*f_cell[1] + coef_fv[2]* f_cell[0] + coef_fv[1]*f_cell[N-1] + coef_fv[0]* f_cell[N-2] 
        df_pred_faceFV[1]   = coef_fv[3]*f_cell[2] + coef_fv[2]* f_cell[1] + coef_fv[1]*f_cell[0]   + coef_fv[0]* f_cell[N-1]
        df_pred_faceFV[N-1] = df_pred_faceFV[1]  
        df_pred_faceFV[N]   = df_pred_faceFV[0]

        f_pred_faceFV[0]   = int_fv[3]*f_cell[1] + int_fv[2]* f_cell[0] + int_fv[1]*f_cell[N-1] + int_fv[0]* f_cell[N-2] 
        f_pred_faceFV[1]   = int_fv[3]*f_cell[2] + int_fv[2]* f_cell[1] + int_fv[1]*f_cell[0]   + int_fv[0]* f_cell[N-1]
        
        f_pred_faceFV[N-1] = int_fv[3]*f_cell[0] + int_fv[2]* f_cell[N-1] + int_fv[1]*f_cell[N-2]   + int_fv[0]* f_cell[N-3]
        f_pred_faceFV[N]   = f_pred_faceFV[0]

    elif (order==6): 
        df_pred_face[0]   = coef_fd[5]*f_vals[2] + coef_fd[4]*f_vals[1] +  coef_fd[3]*f_vals[0]  + coef_fd[2]* f_vals[N-1] + coef_fd[1]*f_vals[N-2] + coef_fd[0]* f_vals[N-3] 
        df_pred_face[1]   = coef_fd[5]*f_vals[3] + coef_fd[4]*f_vals[2] +  coef_fd[3]*f_vals[1]  + coef_fd[2]* f_vals[0]   + coef_fd[1]*f_vals[N-1] + coef_fd[0]* f_vals[N-2] 
        df_pred_face[2]   = coef_fd[5]*f_vals[4] + coef_fd[4]*f_vals[3] +  coef_fd[3]*f_vals[2]  + coef_fd[2]* f_vals[1]   + coef_fd[1]*f_vals[0]   + coef_fd[0]* f_vals[N-1] 
        
        df_pred_face[N-2] = df_pred_face[2]  
        df_pred_face[N-1] = df_pred_face[1]
        df_pred_face[N]   = df_pred_face[0]

        f_pred_face[0]   = int_fd[5]*f_vals[2] + int_fd[4]*f_vals[1] +  int_fd[3]*f_vals[0]  + int_fd[2]* f_vals[N-1] + int_fd[1]*f_vals[N-2] + int_fd[0]* f_vals[N-3] 
        f_pred_face[1]   = int_fd[5]*f_vals[3] + int_fd[4]*f_vals[2] +  int_fd[3]*f_vals[1]  + int_fd[2]* f_vals[0]   + int_fd[1]*f_vals[N-1] + int_fd[0]* f_vals[N-2] 
        f_pred_face[2]   = int_fd[5]*f_vals[4] + int_fd[4]*f_vals[3] +  int_fd[3]*f_vals[2]  + int_fd[2]* f_vals[1]   + int_fd[1]*f_vals[0]   + int_fd[0]* f_vals[N-1] 
       
        f_pred_face[N-2]   = int_fd[5]*f_vals[0] + int_fd[4]*f_vals[N-1] +  int_fd[3]*f_vals[N-2]  + int_fd[2]* f_vals[N-3]   + int_fd[1]*f_vals[N-4]   + int_fd[0]* f_vals[N-5]        
        f_pred_face[N-1]   = int_fd[5]*f_vals[1] + int_fd[4]*f_vals[0]   +  int_fd[3]*f_vals[N-1]  + int_fd[2]* f_vals[N-2]   + int_fd[1]*f_vals[N-3]   + int_fd[0]* f_vals[N-4]        
        f_pred_face[N]     = f_pred_face[0]

        df_pred_faceFV[0]   = coef_fv[5]*f_cell[2] + coef_fv[4]*f_cell[1] +  coef_fv[3]*f_cell[0]  + coef_fv[2]* f_cell[N-1] + coef_fv[1]*f_cell[N-2] + coef_fv[0]* f_cell[N-3] 
        df_pred_faceFV[1]   = coef_fv[5]*f_cell[3] + coef_fv[4]*f_cell[2] +  coef_fv[3]*f_cell[1]  + coef_fv[2]* f_cell[0]   + coef_fv[1]*f_cell[N-1] + coef_fv[0]* f_cell[N-2] 
        df_pred_faceFV[2]   = coef_fv[5]*f_cell[4] + coef_fv[4]*f_cell[3] +  coef_fv[3]*f_cell[2]  + coef_fv[2]* f_cell[1]   + coef_fv[1]*f_cell[0]   + coef_fv[0]* f_cell[N-1] 
        df_pred_faceFV[N-2] = df_pred_faceFV[2]  
        df_pred_faceFV[N-1] = df_pred_faceFV[1]  
        df_pred_faceFV[N]   = df_pred_faceFV[0]

        f_pred_faceFV[0]   = int_fv[5]*f_cell[2] + int_fv[4]*f_cell[1] +  int_fv[3]*f_cell[0]  + int_fv[2]* f_cell[N-1] + int_fv[1]*f_cell[N-2] + int_fv[0]* f_cell[N-3] 
        f_pred_faceFV[1]   = int_fv[5]*f_cell[3] + int_fv[4]*f_cell[2] +  int_fv[3]*f_cell[1]  + int_fv[2]* f_cell[0]   + int_fv[1]*f_cell[N-1] + int_fv[0]* f_cell[N-2] 
        f_pred_faceFV[2]   = int_fv[5]*f_cell[4] + int_fv[4]*f_cell[3] +  int_fv[3]*f_cell[2]  + int_fv[2]* f_cell[1]   + int_fv[1]*f_cell[0]   + int_fv[0]* f_cell[N-1] 
        
        f_pred_faceFV[N-2]   = int_fv[5]*f_cell[0] + int_fv[4]*f_cell[N-1] +  int_fv[3]*f_cell[N-2]  + int_fv[2]* f_cell[N-3]   + int_fv[1]*f_cell[N-4]   + int_fv[0]* f_cell[N-5] 
        f_pred_faceFV[N-1]   = int_fv[5]*f_cell[1] + int_fv[4]*f_cell[0]   +  int_fv[3]*f_cell[N-1]  + int_fv[2]* f_cell[N-2]   + int_fv[1]*f_cell[N-3]   + int_fv[0]* f_cell[N-4] 
        f_pred_faceFV[N]   = f_pred_faceFV[0]

    else:
        raise ValueError("Unsupported order. Only orders 2, 4, and 6 are allowed.")    

    df_pred_face    = df_pred_face/dx
    df_pred_faceFV  = df_pred_faceFV/dx
    
    ########

    
    # Central derivatives at cell-centred   -----------------------
    # derivative at df_cell[i] = 0.5*(f[i+1] -f[i-1]) using 2nd order

    if (order==2):
        for i in range(istart, iend + 1):  
            df_pred[i] = coef_cc[1] *f_vals[i+1] +coef_cc[0]* f_vals[i-1]
        df_pred[0]   = coef_cc[1]*f_vals[1] + coef_cc[0]* f_vals[N-1]
        df_pred[N-1] = coef_cc[1]*f_vals[0] + coef_cc[0]* f_vals[N-2]
    elif (order==4):  
        for i in range(istart, iend + 1):  
            df_pred[i] = coef_cc[3] *f_vals[i+2] + coef_cc[2] *f_vals[i+1] +coef_cc[1] *f_vals[i-1] +coef_cc[0]* f_vals[i-2] 
        df_pred[0]   = coef_cc[3]*f_vals[2] + coef_cc[2]* f_vals[1]   +  coef_cc[1]*f_vals[N-1] + coef_cc[0]* f_vals[N-2] 
        df_pred[1]   = coef_cc[3]*f_vals[3] + coef_cc[2]* f_vals[2]   +  coef_cc[1]*f_vals[0]   + coef_cc[0]* f_vals[N-1] 
        df_pred[N-1] = coef_cc[3]*f_vals[1] + coef_cc[2]* f_vals[0]   +  coef_cc[1]*f_vals[N-2] + coef_cc[0]* f_vals[N-3] 
        df_pred[N-2] = coef_cc[3]*f_vals[0] + coef_cc[2]* f_vals[N-1] +  coef_cc[1]*f_vals[N-3] + coef_cc[0]* f_vals[N-4] 
    elif (order==6):  
        for i in range(istart, iend + 1):  
            df_pred[i] = coef_cc[5] *f_vals[i+3] + coef_cc[4] *f_vals[i+2] + coef_cc[3] *f_vals[i+1] + coef_cc[2] *f_vals[i-1] +coef_cc[1] *f_vals[i-2] +coef_cc[0]* f_vals[i-3] 
        jj = 0
        kk = N 
        df_pred[jj] = coef_cc[5] *f_vals[jj+3] + coef_cc[4] *f_vals[jj+2] + coef_cc[3] *f_vals[jj+1] + coef_cc[2] *f_vals[kk-1] +coef_cc[1] *f_vals[kk-2] +coef_cc[0]* f_vals[kk-3] 
        jj = 1
        df_pred[jj] = coef_cc[5] *f_vals[jj+3] + coef_cc[4] *f_vals[jj+2] + coef_cc[3] *f_vals[jj+1] + coef_cc[2] *f_vals[jj-1] +coef_cc[1] *f_vals[kk-1] +coef_cc[0]* f_vals[kk-2] 
        jj = 2
        df_pred[jj] = coef_cc[5] *f_vals[jj+3] + coef_cc[4] *f_vals[jj+2] + coef_cc[3] *f_vals[jj+1] + coef_cc[2] *f_vals[jj-1] +coef_cc[1] *f_vals[jj-2] +coef_cc[0]* f_vals[kk-1] 

        jj = N-1
        kk = 0
        df_pred[jj] = coef_cc[5] *f_vals[kk+2] + coef_cc[4] *f_vals[kk+1] + coef_cc[3] *f_vals[kk] + coef_cc[2] *f_vals[jj-1] +coef_cc[1] *f_vals[jj-2] +coef_cc[0]* f_vals[jj-3] 
        jj = N-2
        df_pred[jj] = coef_cc[5] *f_vals[kk+1] + coef_cc[4] *f_vals[kk] + coef_cc[3] *f_vals[jj+1] + coef_cc[2] *f_vals[jj-1] +coef_cc[1] *f_vals[jj-2] +coef_cc[0]* f_vals[jj-3] 
        jj = N-3
        df_pred[jj] = coef_cc[5] *f_vals[kk]   + coef_cc[4] *f_vals[jj+2] + coef_cc[3] *f_vals[jj+1] + coef_cc[2] *f_vals[jj-1] +coef_cc[1] *f_vals[jj-2] +coef_cc[0]* f_vals[jj-3] 

    else:
        print(" order 2/4/6 ")
    df_pred = df_pred/dx

    # ------------ ERRORS DERIVATIVE
    print(" --   Derivative Errors -- ")

    # L2 error cell derivative
    abs_error = np.abs(df_pred - df_exact,    dtype=np.float64)
    l2_error = np.sqrt(np.mean(abs_error*abs_error, dtype=np.float64), dtype=np.float64)
    errors.append(l2_error)
    print(f"Resolution: {N}, L2 Error (cell-centred): {l2_error:.5e}")

    # L2 error face derivative FD
    abs_errorfc = np.abs(df_pred_face - df_exact_face,    dtype=np.float64)
    l2_errorfc = np.sqrt(np.mean(abs_errorfc*abs_errorfc, dtype=np.float64), dtype=np.float64)
    errors_fc.append(l2_errorfc)
    print(f"Resolution: {N}, L2 Error (face-centred)  FD: {l2_errorfc:.5e}")

    # L2 error face derivative FV
    abs_errorfv = np.abs(df_pred_faceFV - df_exact_face,    dtype=np.float64)
    l2_errorfv = np.sqrt(np.mean(abs_errorfv*abs_errorfv, dtype=np.float64), dtype=np.float64)
    errors_fv.append(l2_errorfv)
    print(f"Resolution: {N}, L2 Error (face-centred)  FV: {l2_errorfv:.5e}")
    # ------------ ERRORS INTERPOLATION
    print(" --   Interpolation Errors -- ")
    # L2 error FD
    abs_error1 = np.abs(f_pred_face - f_exact_face,    dtype=np.float64)
    l2_error1 = np.sqrt(np.mean(abs_error1*abs_error1, dtype=np.float64), dtype=np.float64)
    errors_int.append(l2_error1)
    print(f"Resolution: {N}, L2 Error (face-centred)  FD: {l2_error1:.5e}")

    # L2 error FV
    abs_error2 = np.abs(f_pred_faceFV - f_exact_face,    dtype=np.float64)
    l2_error2 = np.sqrt(np.mean(abs_error2*abs_error2, dtype=np.float64), dtype=np.float64)
    errors_intfv.append(l2_error2)
    print(f"Resolution: {N}, L2 Error (face-centred)  FV: {l2_error2:.5e}")


# ---------------------------------------------------------------------------------------


# Plot df

# plt.figure(figsize=(8, 6))
# plt.plot(x_cell,df_pred, marker = 'o',label='pred')
# plt.plot(x_cell,df_exact,label='exact')
# plt.title('Cell Centered Derivative')
# plt.legend()
# plt.show()


# uncoment to see Derivative

# plt.figure(figsize=(8, 6))
# plt.plot(x_face,df_pred_face, marker = 'o',label='pred FD')
# plt.plot(x_face,df_pred_faceFV, marker = '+',label='pred FV')
# plt.plot(x_face,df_exact_face,label='exact')
# plt.title('Face Centered Derivative')
# plt.legend()
# plt.show()


# plt.figure(figsize=(8, 6))
# plt.plot(x_cell,f_cell, marker = 'o',label='FV')
# plt.plot(x_cell,f_vals,label='exact (FD)')
# plt.title('Evaluated Function FD or FV ')
# plt.legend()
# plt.show()

plt.figure(figsize=(8, 6))
plt.plot(x_face,f_pred_face, marker = 'o',label='pred FD')
plt.plot(x_face,f_pred_faceFV, marker = '+',label='pred FV')
plt.plot(x_face,f_exact_face,label='exact')
plt.title('Evaluated Function at cell faces FD or FV ')
plt.legend()
plt.show()




# Plot convergence

plt.figure(figsize=(8, 6))
e0    = errors[0]
e0_fc = errors_fc[0]
e0_fv = errors_fv[0]

plt.loglog(resolutions, errors/e0, 'x', color='green',label='L2 Error Cell-centred')

plt.loglog(resolutions, errors_fc/e0_fc, 'o', color='orange', label='L2 Error Face-centred (FD)')

plt.loglog(resolutions, errors_fv/e0_fv, '+', color='blue', label='L2 Error Face-centred (FV)')



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
plt.title('Convergence of First Derivative')
plt.legend()

# Optional: fit and plot convergence rate
p = np.polyfit(np.log(resolutions), np.log(errors), 1)
order = -p[0]      
plt.text(resolutions[0], errors[0], f"Order ≈ {order:.2f}", fontsize=12)

plt.show()


plt.figure(figsize=(8, 6))
ei0    = errors_int[0]
ei0_fv = errors_intfv[0]

plt.loglog(resolutions, errors_int/ei0, 'o', color='orange',label='L2 Error using FD')
plt.loglog(resolutions, errors_intfv/ei0_fv, '+', color='blue', label='L2 Error using FV')

plt.loglog(resolutions, ideal_2nd, linestyle='--', color='gray', label='2nd-order')
plt.loglog(resolutions, ideal_4th, linestyle='--', color='gray', label='4th-order')
plt.loglog(resolutions, ideal_6th, linestyle='--', color='gray', label='6th-order')


plt.grid(True, which='both', ls='--')
plt.xlabel('Resolution (N)')
plt.ylabel('Normalised L2 Error')
plt.title('Convergence of Interpolation Face')
plt.legend()

# Optional: fit and plot convergence rate
p = np.polyfit(np.log(resolutions), np.log(errors_int), 1)
order = -p[0]      
plt.text(resolutions[0], errors_int[0], f"Order ≈ {order:.2f}", fontsize=12)

plt.show()


# Example results

#  --   Interpolation Errors --  Order 2
# Resolution: 128, L2 Error (face-centred)  FD: 2.12140e-04
# Resolution: 128, L2 Error (face-centred)  FV: 2.82834e-04

#  --   Interpolation Errors --  Order 4
#Resolution: 128, L2 Error (face-centred)  FD: 1.21367e-07
#Resolution: 128, L2 Error (face-centred)  FV: 2.87664e-07

#  --   Derivative Errors --   Order 2
# Resolution: 128, L2 Error (cell-centred): 1.78403e-03
# Resolution: 128, L2 Error (face-centred)  FD: 4.47772e-04
# Resolution: 128, L2 Error (face-centred)  FV: 8.95500e-04

# --   Derivative Errors -- 
#Resolution: 128, L2 Error (cell-centred): 4.43764e-10
#Resolution: 128, L2 Error (face-centred)  FD: 4.35101e-11
#Resolution: 128, L2 Error (face-centred)  FV: 1.11390e-10
# --   Interpolation Errors -- 
#Resolution: 128, L2 Error (face-centred)  FD: 4.80973e-11
#Resolution: 128, L2 Error (face-centred)  FV: 7.03530e-11



