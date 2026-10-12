import os
from dotenv import load_dotenv
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import windows, csd

def postprocess(Uinput, parameters):
    U = np.loadtxt(Uinput, skiprows=parameters.uskip, delimiter=',')
    
    #sampling frequency
    fsample = 1 / (U[1,0] - U[0,0])
    
    #domain flow-through time
    #NOTE: different normalisation used for fft. 
    dft = U[:,0] / (parameters.Lx / parameters.u_ref) 
    dft_fft = U[:,0] / (parameters.L_ref / parameters.u_ref) 
    
    #yvelocity computation from ymom and rho
    #vel normalised by freestream velocity
    Ymom = U[:,2]
    density = U[:,3]
    yvel = Ymom / density
    yvel /= parameters.u_ref

    plt.plot(dft, yvel, label='u_y', color='green', linewidth=2)
    plt.ylabel('v/Uinf')
    plt.xlabel('t*U/Lx')
    plt.show()

    return yvel, dft_fft, fsample

def fft(Upp, t_star, fs, parameters):
    nwindows = 1
    nff = None

    #VERIFY both non-normalisations produce same results.
    fs = fs * parameters.L_ref / parameters.u_ref
    #fs = 1/ (t_star[2] - t_star[1])

    f, Puu = csd(Upp, Upp, fs, window="hann", nperseg=int(len(Upp) / nwindows), nfft=nff)
    plt.plot(f, Puu, label='Puu', color='red', linewidth=2.5)
    plt.yscale('log')
    plt.ylim(10**-8, 5)
    plt.xlim((0, 3))
    plt.ylabel('v/Uinf (power spectral density)')
    plt.xlabel('f*L/u')
    plt.show()

class Parameters:
    def __init__(self, run):
        self.uskip = 4

        if run == "M6":
            self.u_ref = 208.3500981
            self.L_ref = 0.0002256
            self.cL = 13.184
        elif run == "M5":
            self.u_ref = 173.6250818
            self.L_ref = 0.0002708
            self.cL = 13.662
        else:
            raise ValueError(f"Unknown run configuration: {run}")
        
        self.Lx = self.cL * self.L_ref / 2 

if __name__ == "__main__":
    # 1. Read configurations from environment variables (with safe repo defaults)
    load_dotenv()
    RUN_CASE = os.getenv("SIM_RUN", None)
    DATA_ROOT = Path(os.getenv("SIM_DATA_DIR", None))
    SIGNAL_NAME = os.getenv("SIM_SIGNAL_NAME", None)
    CAVS = os.getenv("SIM_CAV", None)

    # 2. Comma-separated cavity list or single default
    params = Parameters(RUN_CASE)

    # 3. Robust cross-platform path construction
    input_file = DATA_ROOT / CAVS / f"{SIGNAL_NAME}.log"
    
    if not input_file.exists():
        print(f"File not found: {input_file}")

    U_norm, t_star, fsample = postprocess(input_file, params)
    fft(U_norm, t_star, fsample, params)