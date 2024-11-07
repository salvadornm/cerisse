import yt
import numpy as np
import matplotlib.pyplot as plt
import glob
import os

Lx = 0.08214
Ly = 0.01369
Lz = 0.02738
xc = 0.5*Lx
yc = 0.5*Lx
zc = 0.5*Lz


#####################
def compute_wall_units(y, u, tau_w, rho, mu):
  u_tau = np.sqrt(tau_w/rho)
  nu = mu/rho
  # compute y+ u+
  y_plus = y *u_tau/nu
  u_plus = u / u_tau

  return y_plus, u_plus
##################
def ExtractLineoutY(datasets, xc,zc):
  y_out = []
  v_out = []
  
  for d in range(len(datasets)):
    ds = yt.load(datasets[d])

    axis = 1 # take a line cut along y axis  (0: x  2:z)
    
    zp = zc
    xp = xc
 
    lineout = ds.ortho_ray(axis,(zp, xp))
    
    srt = np.argsort(lineout["y"])
    y = np.array(lineout["y"][srt])
    rhov = np.array(lineout["Xmom"][srt])
    rho  =  np.array(lineout["Density"][srt])
    v = rhov/rho 

    y_out.append(y)
    v_out.append(v)

  return y_out, v_out
################################################

print('opening data')

# open data
list_of_files = glob.glob('./plot/*')
latest_file = max(list_of_files, key=os.path.getctime)
print(" LAST FILE=",latest_file)
datasets = [latest_file]

y, v  = ExtractLineoutY(datasets, xc, zc)


for (ys, vs) in zip(y, v): 
  
  # Exact  solution
  #plt.plot(ve,ys,'o',label='Analytical',mfc='none',markersize=3)

  # cerisse
  plt.plot(vs, ys, label='cerisse')
  plt.ylabel(' y ')
  plt.xlabel(' v ', horizontalalignment='center')#####################
def ExactSolution(y):
  yoh = y/L
  u = 1.5*Vbulk*(4*(1 - yoh)*yoh)  
  return u

  plt.legend()
  plt.show()

# Save the line plot into a file
#plot.save()


print(' .... DONE')






print(' .... DONE')
