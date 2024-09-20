import yt
import numpy as np
import matplotlib.pyplot as plt
import glob
import os

L = 0.3112
R = L/20  # vortex size
v0 = 35.0
beta = 0.04  # vortex strength

Gam  = beta*v0*R*np.sqrt(np.exp(1.0))
xc = 0.5*L
yc = 0.5*L
rho0 = 1.17170407

print(" Vortex parameters ..")

print(" Gam ",Gam," beta=",beta)


#####################
def ExactSolution(x):
  u = -v0*beta*(x-yc)/R*np.exp(-0.5*((x-yc)/R)**2)
  v = 2*Gam/(R*R)*np.exp(-((x-xc)/R)**2)*(x-xc)  

  return u,v
##################
def VortexLineoutX(datasets, xc, yc):
  x_out = []
  v_out = []
  rho_out = []

  for d in range(len(datasets)):
    ds = yt.load(datasets[d])

    xaxis = 0  # take a line cut along x axis
    lineout = ds.ortho_ray(xaxis, (yc, 0))
    srt = np.argsort(lineout["x"])
    x = np.array(lineout["x"][srt])
    v = np.array(lineout["Ymom"][srt])
    x_out.append(x)
    v_out.append(v)
    rho = np.array(lineout["Density"][srt])
    rho_out.append(rho)

  return x_out, v_out, rho_out
################################################

print('opening data')

# open data
list_of_files = glob.glob('./plot/*') 
latest_file = max(list_of_files, key=os.path.getctime)
print(" LAST FILE=",latest_file)
datasets = [latest_file]

#print("extracting lines at xc,yc=",xc,yc)

x, v, rho = VortexLineoutX(datasets, xc, yc)

for (xs, vs, rhos) in zip(x, v, rho): 
  plt.plot(xs, vs/rhos, label='cerisse')
  plt.ylabel("vy")
  plt.xlabel(' x ', horizontalalignment='center')

  ue,ve  = ExactSolution(xs)
  
  # exact solution
  plt.plot(xs,ve,'o',label='exact')
  plt.legend()
  plt.show()

# Save the line plot into a file
#plot.save()


print(' .... DONE')






print(' .... DONE')
