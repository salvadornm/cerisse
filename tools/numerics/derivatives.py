import numpy as np

order = 4  #<============  2/4/6
HALO = int(order/2)

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

stencil = 2*HALO
print(" ----------  ")
print(" Order=", order)
print(" stencil=",stencil, " halo ",HALO)
print(" ----------  ")

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

iv = [1, 0, 0]
jv = [0, 1, 0]
kv = [0, 0 ,1]


# CENTRAL DERIVATIVES CELL_CENTRED  df = 0.5*(f(i+1) - f(i-1))
# derivative of f in direction idir at cell i,j,k, asuming mesh spacing of 1
# idir  : 0,1,2
# i j k : cell
# f: numpy 3D array
def dfdx_cc(i,j,k,idir,f):
    ii = iv[idir] 
    jj = jv[idir] 
    kk = kv[idir]
    i1 = i - HALO*ii   # i.e:  idir=0  order = 4 (HALO=2) i1 = i-2
    j1 = j - HALO*jj
    k1 = k - HALO*kk

    df = 0.0
    for n in range(0, HALO):   # i.e: order 4 looping n=[0,2]  n=0,1
        # print(" n=",n, " coef[n]=",coef_cc[n]," coef[n+HALO]=", coef_cc[n+HALO]) 
        # print(" i -1 ",i1 + n*ii, " i+1",i  + (n+1)*ii)
        df += coef_cc[n]*f[i1 + n*ii, j1 + n*jj, k1  +n*kk]                   #  i.e.:  coef[0]*f[i-2](n=0) coef[1]*f[i-1](n=1)
        df += coef_cc[n+HALO]*f[i  + (n+1)*ii, j  + (n+1)*jj, k   +(n+1)*kk]  #  i.e.:  coef[2]*f[i+1](n=0) coef[3]*f[i+2](n=1)
    return df


# Interpolation  f[i+1/2]= 0.5*(f(i) + f(i+1)) 
# inteporlation of f in  face i+1/2  FINITE DIFFERENCE style 
# idir  : 0,1,2
# nc    : 0,1,2,3,4
# i j k : cell
# f: numpy 3D array
def interp_fd(i,j,k,idir,nc,f):
    ii = iv[idir] 
    jj = jv[idir] 
    kk = kv[idir]
    i1 = i - HALO*ii   # i.e:  idir=0  order = 4 (HALO=2) i1 = i-2
    j1 = j - HALO*jj
    k1 = k - HALO*kk

    df = 0.0
    for n in range(0, stencil):   # i.e: order 4 looping n=[0,2]  n=0,1
        df += int_fd[n]*f[i1 + n*ii, j1 + n*jj, k1  +n*kk,nc]               
    return df

# Interpolation  f[i+1/2]= 0.5*(f(i) + f(i+1)) 
# inteporlation of f in  face i+1/2  FINITE VOLUME style 
# idir  : 0,1,2
# nc    : 0,1,2,3,4
# i j k : cell
# f: numpy 3D array
def interp_fv(i,j,k,idir,nc,f):
    ii = iv[idir] 
    jj = jv[idir] 
    kk = kv[idir]
    i1 = i - HALO*ii   # i.e:  idir=0  order = 4 (HALO=2) i1 = i-2
    j1 = j - HALO*jj
    k1 = k - HALO*kk

    df = 0.0
    for n in range(0, stencil):   # i.e: order 4 looping n=[0,2]  n=0,1
        df += int_fv[n]*f[i1 + n*ii, j1 + n*jj, k1  +n*kk,nc]               
    return df
# Interpolation  phi[i+1/2]= 0.5*(phi(i) + phi(i+1)) 
# interpolation of phi in  face i+1/2  FINITE VOLUME style 
def interp_fv_phi(i,j,k,idir,phi):
    ii = iv[idir] 
    jj = jv[idir] 
    kk = kv[idir]
    i1 = i - HALO*ii   # i.e:  idir=0  order = 4 (HALO=2) i1 = i-2
    j1 = j - HALO*jj
    k1 = k - HALO*kk

    df = 0.0
    for n in range(0, stencil):   # i.e: order 4 looping n=[0,2]  n=0,1
        df += int_fv[n]*phi[i1 + n*ii, j1 + n*jj, k1  +n*kk]               
    return df

# Derivative  df[i+1/2]= f(i+1) - f(i)
# deivative at face i+1/2  FINITE VOLUME style 
# idir  : 0,1,2
# i j k : cell
# f: numpy 3D array
def dfdx_fv(i,j,k,idir,f):
    ii = iv[idir] 
    jj = jv[idir] 
    kk = kv[idir]
    i1 = i - HALO*ii  
    j1 = j - HALO*jj
    k1 = k - HALO*kk

    df = 0.0
    for n in range(0, stencil):   
        df += coef_fv[n]*f[i1 + n*ii, j1 + n*jj, k1  +n*kk]               
    return df

# Derivative  df[i+1/2]= f(j+1) - f(j)
def dfcross_fv(i,j,k,idir,idir2,f):
    ii = iv[idir] 
    jj = jv[idir] 
    kk = kv[idir]

    i1 = i - HALO*ii  
    j1 = j - HALO*jj
    k1 = k - HALO*kk

    # derivatives cell centred
    df = 0
    # loop over cells isten
    for n in range(0, stencil): 
        d1y =   dfdx_cc(i1 + n*ii,j1+ n*jj,k1+ n*kk,idir2,f)
        df += int_fv[n]*d1y            
    return df    


def test(i,j,k,idir,f):

    aux = 67
    print(" i j k ",i,j,k)
    ii = iv[idir] 
    jj = jv[idir] 
    kk = kv[idir]
    print(" ii jj kk ",ii,jj,kk)
    print(" f = ", f[i,j,k])

    return aux
