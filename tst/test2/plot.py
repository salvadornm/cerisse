import yt
import numpy as np
import matplotlib.pyplot as plt
#-------------------------------------------
# file to open data and plot 1D results
# tst/tst1
#-----------------------------------------

print('opening data')

for fname, color in zip(["plot/plt00595", "plot/plt00953"], ["C0", "C1",]):
  # load data
  ds = yt.load(fname)

  xaxis = 0  # take a line cut along the x axis
  lineout = ds.ortho_ray(xaxis, (0, 0))   # cutting through the xaxis
  srt = np.argsort(lineout["index", "x"]) # sort according to x position

  # plot solutions
  for i, name in zip([1, 2, 3], ["temperature", "pressure", "rho_AR"]):
    plt.subplot(3, 1, i)
    plt.plot(np.array(lineout["index", "x"][srt]),np.array(lineout["boxlib", name][srt]), label=rf"{float(ds.current_time)*1e6:.0f} $\mu$s", color=color)
    plt.ylabel(name)
    plt.xlim([float(ds.domain_left_edge[0]), float(ds.domain_right_edge[0])])
ref_data = np.loadtxt("ref.dat", skiprows=1) 
plt.subplot(3, 1, 1)
plt.plot(ref_data[:,0], ref_data[:,1], 'o', color="C0")
plt.plot(ref_data[:,0], ref_data[:,4], 'o', color="C1")
plt.subplot(3, 1, 2)
plt.plot(ref_data[:,0], ref_data[:,2], 'o', color="C0")
plt.plot(ref_data[:,0], ref_data[:,5], 'o', color="C1")
plt.legend()
plt.subplot(3, 1, 3)
plt.plot(ref_data[:,0], ref_data[:,3] * 0.88586543, 'o', color="C0")
plt.plot(ref_data[:,0], ref_data[:,6] * 0.88586543, 'o', color="C1")

plt.xlabel("x [m]")
plt.tight_layout()

### plt.legend()
plt.show()

#plt.savefig("test2_ref.png", facecolor="white")

print(' .... DONE')
