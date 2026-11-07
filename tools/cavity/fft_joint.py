import numpy as np
import os
from dataclasses import dataclass
import matplotlib.pyplot as plt
import scipy as sp
from scipy.signal import windows, csd

def postprocess(Uinput, Uout, parameters, i, label_name):
    U = np.loadtxt(Uinput, skiprows=parameters.uskip, delimiter=',')
    #sampling frequency
    fsample = 1 / (U[1,0] - U[0,0])
    #domain flow-through time
    #dft = U[:,0]
    dft = U[:,0] / (parameters.Lx[i] / parameters.u_ref)
    #dft = U[:,0] / (parameters.L_ref / parameters.u_ref)
    #yvelocity
    Ymom = U[:,2]
    density = U[:,3]
    yvel = Ymom / density
    yvel /= parameters.u_ref

    # Combine the time and velocity arrays into two columns
    output_data = np.column_stack((dft, yvel))

    # Write to the specified output file
    np.savetxt(Uout, output_data, delimiter=',', header='t*U/Lx,v/Uinf', comments='')
    
    # Target Figure 1 for time-series and use dynamic colors/labels
    plt.figure(1)
    plt.plot(dft, yvel, label=label_name, linewidth=2)
    plt.ylabel('v/Uinf')
    plt.xlabel('t*U/Lx')

    return yvel, fsample

def fft(Upp, fs, parameters, label_name):
    nwindows = 1
    nff = None
    fs = fs * parameters.L_ref / parameters.u_ref
    f, Puu = csd(Upp, Upp, fs, window='hann', nperseg=len(Upp)/nwindows, nfft=nff)
    
    # Target Figure 2 for FFT and use dynamic colors/labels
    plt.figure(2)
    plt.plot(f, Puu, label=label_name, linewidth=2.5)
    plt.yscale('log')
    plt.ylim(10**-8, 10**0)
    plt.xlim((0, 3))
    plt.ylabel('v/Uinf (power spectral density)')
    plt.xlabel('f*L/u')
    return

class Parameters:
    def __init__(self):
        self.uskip = 4
        self.u_ref = 208.3500981
        self.L_ref = 0.0002256
        #self.cL = [13.184, 16.324, 19.465]
        self.cL = [13.184,13.184,13.184]
        self.Lx = [x * self.L_ref / 2 for x in self.cL]

params = Parameters()

#cavs = ["cav4_theta2", "cavity_thetaB", "cavity_thetaBB"]
cavs = ["sigma01", "sigma001", "sigma005"]

# Process all files and add lines to the figures
for i, cav in enumerate(cavs):
    base_path = f"./../../exm/ebm/reflection/{cav}"
    name = "time_probe_react"
    U_in = f"{base_path}/{name}.log"
    U_ou = f"{base_path}/{name}_out.log"
    
    # Pass the cavity name (cav) as the label
    U_norm, fsample = postprocess(U_in, U_ou, params, i, cav)
    fft(U_norm, fsample, params, cav)

# Finalize and show the plots after the loop is done
plt.figure(1)
plt.legend()
plt.title("Time History")

plt.figure(2)
plt.legend()
plt.title("FFT - Power Spectral Density")

plt.show()
