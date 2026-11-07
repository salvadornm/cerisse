import numpy as np
import os
from dataclasses import dataclass
import matplotlib.pyplot as plt
import scipy as sp
from scipy.signal import windows, csd

def postprocess(Uinput,Uout,parameters,i):
    U = np.loadtxt(Uinput,skiprows=parameters.uskip,delimiter=',')
    #sampling frequency
    fsample = 1 / (U[1,0] - U[0,0])
    #domain flow-through time
    dft = U[:,0] / (parameters.Lx[i] / parameters.u_ref)
    #yvelocity
    Ymom = U[:,2]
    density = U[:,3]
    yvel = Ymom / density
    yvel /= parameters.u_ref

    # Combine the time and velocity arrays into two columns
    output_data = np.column_stack((dft, yvel))

# Write to the specified output file
    np.savetxt(Uout, output_data, delimiter=',', header='t*U/Lx,v/Uinf', comments='')
    plt.plot(dft, yvel, label='u_y', color='green', linewidth=2)
    plt.ylabel('v/Uinf')
    plt.xlabel('t*U/Lx')
    plt.show()

    return yvel,fsample

def fft(Upp,fs,parameters):
    nwindows = 1
    nff = None
    fs = fs* parameters.L_ref / parameters.u_ref
    f, Puu = csd(Upp,Upp,fs,window='hann',nperseg= len(Upp)/nwindows,nfft=nff)
    #Puu = Puu * (parameters.u_ref / parameters.L_ref)
    #f = f*parameters.L_ref / parameters.u_ref
    plt.plot(f, Puu, label='Puu', color='red', linewidth=2.5)
    plt.yscale('log')
    plt.ylim(10**-8, 10**0)
    plt.xlim((0,3))
    plt.ylabel('v/Uinf (power spectral density)')
    plt.xlabel('f*L/u')
    plt.show()
    return

#L/D = 2
#cavity theta Lx = 13.184D 
#cavity theta B Lx = 16.324D
#cavity theta BB Lx = 19.465D

class Parameters:
    def __init__(self):
        self.uskip = 4
        self.u_ref = 208.3500981
        self.L_ref = 0.0002256
        self.cL = [13.184, 16.324, 19.465]
        
        # Calculate Lx using a simple list comprehension
        self.Lx = [x * self.L_ref / 2 for x in self.cL]
params = Parameters()

cavs = ["cav4_theta2","cavity_thetaB","cavity_thetaBB"]

for i,cav in enumerate(cavs):
    base_path = f"./../../exm/ebm/{cav}/cx3"
    name = "tstart_15"
    U_in = f"{base_path}/{name}.log"
    U_ou = f"{base_path}/{name}_out.log"
    U_norm,fsample = postprocess(U_in,U_ou,params,i)
    fft(U_norm,fsample,params)