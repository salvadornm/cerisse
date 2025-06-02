import sympy as sp
import numpy as np
import matplotlib.pyplot as plt

# Symbolic function for exact derivative
x = sp.symbols('x')
f_sym = sp.sin(2 * sp.pi * x)
df_sym = sp.diff(f_sym, x)
f_func = sp.lambdify((x), f_sym, modules=["numpy"])
df_func = sp.lambdify((x), df_sym, modules=["numpy"])

# Resolutions to test (finer will not work)
resolutions = [16, 32, 64, 128, 256, 512]
#resolutions = [100]

errors = []

errors_fc = []

#  derivatives coefficients cell centred
ccw2_0  = -0.5
ccw2_1  =  0.5

ccw4_0  =   1.0 / 12.0 
ccw4_1  =  -2.0 / 3.0 
ccw4_2  =   2.0 / 3.0 
ccw4_3  =  -1.0 / 12.0 

ccw6_0  =  -1.0 / 60.0 
ccw6_1  =   3.0 / 20.0 
ccw6_2  =  -3.0 / 4.0 
ccw6_3  =   3.0 / 4.0; 
ccw6_4  =  -3.0 / 20.0 
ccw6_5  =  +1.0 / 60.0 

#  derivatives coefficients face centred

fcw2_0 =  -1.0 
fcw2_1 =   1.0

coef_fc2 =  [-1, 1]

fcw4_0 =   1.0 / 24.0; 
fcw4_1 =  -9.0 / 8.0; 
fcw4_2 =   9.0 / 8.0; 
fcw4_3 =  -1.0 / 24.0; 

coef_fc4 =  [1/24, -9/8, 9/8, -1/24]

fcw6_0 =  -1.0  / 60.0; 
fcw6_1 =  +9.0  / 60.0; 
fcw6_2 =  -45.0 / 60.0; 
fcw6_3 =  +45.0 / 60.0; 
fcw6_4 =  -9.0  / 60.0; 
fcw6_5 =  +1.0  / 60.0; 

# [Fraction(-3, 640), Fraction(25, 384), Fraction(-75, 64), Fraction(75, 64), Fraction(-25, 384), Fraction(3, 640)]

coef_fc6 = [-3/640, 25/384, -75/64, 75/64, -25/384, 3/640]

#  interpolation face centred
int_fcw2_0 =  0.5; 
int_fcw2_1 =  0.5; 

int_fcw4_0 = -1.0 / 12.0; 
int_fcw4_1 =  7.0 / 12.0; 
int_fcw4_2 =  7.0 / 12.0; 
int_fcw4_3 = -1.0 / 12.0; 

int_fcw6_0 =  1.0  /  60.0; 
int_fcw6_1 =  -8.0 / 60.0; 
int_fcw6_2 =  37.0 / 60.0; 
int_fcw6_3 =  37.0 / 60.0; 
int_fcw6_4 =  -8.0 / 60.0; 
int_fcw6_5 =  1.0  / 60.0; 

# domain size
L = 1.0

order = 6  #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
HALO = int(order/2)
stencil = 2*HALO
print(" stencil=",stencil, " halo ",HALO)
coef_cc = np.zeros(stencil)
coef_fc = np.zeros(stencil)
int_fc = np.zeros(stencil)

if (order==2):
    coef_cc[0] = ccw2_0
    coef_cc[1] = ccw2_1

    coef_fc[0] = fcw2_0
    coef_fc[1] = fcw2_1

    int_fc[0]  = int_fcw2_0
    int_fc[1]  = int_fcw2_1
    
elif (order==4):
    coef_cc[0] = ccw4_0
    coef_cc[1] = ccw4_1
    coef_cc[2] = ccw4_2
    coef_cc[3] = ccw4_3

    coef_fc = coef_fc4

    int_fc[0]  = int_fcw4_0
    int_fc[1]  = int_fcw4_1
    int_fc[2]  = int_fcw4_2
    int_fc[3]  = int_fcw4_3

elif (order==6):
    coef_cc[0] = ccw6_0
    coef_cc[1] = ccw6_1
    coef_cc[2] = ccw6_2
    coef_cc[3] = ccw6_3
    coef_cc[4] = ccw6_4
    coef_cc[5] = ccw6_5

    coef_fc = coef_fc6 

    int_fc[0]  = int_fcw6_0
    int_fc[1]  = int_fcw6_1
    int_fc[2]  = int_fcw6_2
    int_fc[3]  = int_fcw6_3
    int_fc[4]  = int_fcw6_4
    int_fc[5]  = int_fcw6_5
    
else:
   raise ValueError("Unsupported order. Only orders 2, 4, and 6 are allowed.")
    


for N in resolutions:
    dx = L / N

    # cells are centred internal point   x[0]= dx/2  N cells
    x_cell = np.zeros(N)
    for i in range(0, N):  
        x_cell[i] =  i*dx + 0.5*dx

    # cells are centred internal point   x[0]= 0  N+1 faces
    x_face = np.zeros(N+1)
    for i in range(0, N+1):  
        x_face[i] =  i*dx

    # exact derivative
    df_exact = df_func(x_cell)

    # exact derivative in faces
    df_exact_face = df_func(x_face)
    
    # function evaluated at cell centred
    f_vals = f_func(x_cell)

    df_pred      = np.zeros(N)
    df_pred_face = np.zeros(N+1)
        
    istart = HALO        # first inner cell point   1/2/3
    iend   = N-HALO -1   # last  inner cell point   N-2/N-3/N-4

    # Central derivatives at face-centred   -----------------------
    # derivative at df_face[i] = f[i] -f[i-1] 

    #print(" istart = ",istart," iend = ",iend) 
    for i in range(istart, iend + 2):
        for j in range(0, stencil):  # 0,1
            df_pred_face[i] = df_pred_face[i] + coef_fc[j]*f_vals[i-HALO+j]
   
    # BC 
    if (order==2):
        df_pred_face[0]   = coef_fc[1]*f_vals[0] + coef_fc[0]* f_vals[N-1]
        df_pred_face[N]   = df_pred_face[0]
    elif (order==4):  
        df_pred_face[0]   = coef_fc[3]*f_vals[1] + coef_fc[2]* f_vals[0] + coef_fc[1]*f_vals[N-1] + coef_fc[0]* f_vals[N-2] 
        df_pred_face[1]   = coef_fc[3]*f_vals[2] + coef_fc[2]* f_vals[1] + coef_fc[1]*f_vals[0]   + coef_fc[0]* f_vals[N-1]
        df_pred_face[N-1] = df_pred_face[1]  
        df_pred_face[N]   = df_pred_face[0]
    elif (order==6): 
        df_pred_face[0]   = coef_fc[5]*f_vals[2] + coef_fc[4]*f_vals[1] +  coef_fc[3]*f_vals[0]  + coef_fc[2]* f_vals[N-1] + coef_fc[1]*f_vals[N-2] + coef_fc[0]* f_vals[N-3] 
        df_pred_face[1]   = coef_fc[5]*f_vals[3] + coef_fc[4]*f_vals[2] +  coef_fc[3]*f_vals[1]  + coef_fc[2]* f_vals[0]   + coef_fc[1]*f_vals[N-1] + coef_fc[0]* f_vals[N-2] 
        df_pred_face[2]   = coef_fc[5]*f_vals[4] + coef_fc[4]*f_vals[3] +  coef_fc[3]*f_vals[2]  + coef_fc[2]* f_vals[1]   + coef_fc[1]*f_vals[0]   + coef_fc[0]* f_vals[N-1] 
        df_pred_face[N-2] = df_pred_face[2]  
        df_pred_face[N-1] = df_pred_face[1]  
        df_pred_face[N]   = df_pred_face[0]
    else:
        raise ValueError("Unsupported order. Only orders 2, 4, and 6 are allowed.")    

    df_pred_face  = df_pred_face/dx
    ########

    
    # Central derivatives at cell-centred   -----------------------
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


    # i=50

    # print(" xface[i]=",x_face[i])
    # print(" xcell[i-1] xcel[i]=",x_cell[i-1],x_cell[i])
    # print(" f[i-1] f[i]",f_vals[i-1],f_vals[i])
    # df2  = coef_fc2[1]*f_vals[i]   + coef_fc2[0]* f_vals[i-1] 
    # df4  = coef_fc4[3]*f_vals[i+1] + coef_fc4[2]* f_vals[i]   + coef_fc4[1]*f_vals[i-1] + coef_fc4[0]* f_vals[i-2] 
    # df6  = coef_fc6[5]*f_vals[i+2] + coef_fc6[4]* f_vals[i+1] + coef_fc6[3]*f_vals[i]  + coef_fc6[2]* f_vals[i-1] + coef_fc6[1]*f_vals[i-2] + coef_fc6[0]* f_vals[i-3] 
    # print(" dfdx (2)",df2/dx)
    # print(" dfdx (4)",df4/dx)
    # print(" dfdx (6)",df6/dx)

    # print(" dfpred_face[i]=",df_pred_face[i])
    # print(" dfexact_face[i]=",df_exact_face[i])
    # print(" exact = ",2*np.pi*np.cos(2*np.pi*x_face[i]))


    # ---------------------------------------------------------
    # L2 error
    abs_error = np.abs(df_pred - df_exact,    dtype=np.float64)
    l2_error = np.sqrt(np.mean(abs_error**2, dtype=np.float64), dtype=np.float64)
    errors.append(l2_error)
    print(f"Resolution: {N}, L2 Error (cell-centred): {l2_error:.5e}")

    abs_errorfc = np.abs(df_pred_face - df_exact_face,    dtype=np.float64)
    l2_errorfc = np.sqrt(np.mean(abs_errorfc**2, dtype=np.float64), dtype=np.float64)
    errors_fc.append(l2_errorfc)
    print(f"Resolution: {N}, L2 Error (face-centred): {l2_errorfc:.5e}")


# ---------------------------------------------------------------------------------------


# Plot df

# plt.figure(figsize=(8, 6))
# plt.plot(x_cell,df_pred, marker = 'o',label='pred')
# plt.plot(x_cell,df_exact,label='exact')
# plt.title('Cell Centered Derivative')
# plt.legend()
# plt.show()


plt.figure(figsize=(8, 6))
plt.plot(x_face,df_pred_face, marker = 'o',label='pred')
plt.plot(x_face,df_exact_face,label='exact')
plt.title('Face Centered Derivative')
plt.legend()
plt.show()


# Plot convergence

plt.figure(figsize=(8, 6))
e0    = errors[0]
e0_fc = errors_fc[0]

plt.loglog(resolutions, errors/e0, marker='o', label='L2 Error')

plt.loglog(resolutions, errors_fc/e0_fc, marker='o', label='L2 Error Face-centred')


# Add ideal convergence lines
resolutions = np.array(resolutions)  # Convert list to NumPy array
N_ref = resolutions[0]
ideal_2nd = (resolutions / N_ref) ** -2
ideal_4th = (resolutions / N_ref) ** -4

#plt.loglog(resolutions, ideal_2nd, 'k--', label='2nd-order')
#plt.loglog(resolutions, ideal_4th, 'k-.', label='4th-order')

plt.grid(True, which='both', ls='--')
plt.xlabel('Resolution (N)')
plt.ylabel('L2 Error')
plt.title('Convergence of Cell Centered Derivative')
plt.legend()

# Optional: fit and plot convergence rate
p = np.polyfit(np.log(resolutions), np.log(errors), 1)
order = -p[0]      
plt.text(resolutions[0], errors[0], f"Order ≈ {order:.2f}", fontsize=12)


plt.show()
