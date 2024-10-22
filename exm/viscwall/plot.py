import yt
import numpy as np
import matplotlib.pyplot as plt
import glob
import os

L = 0.016
beta = 0.04  # vortex strength

xc = 0.5*L
yc = 0.5*L
T0 = 300.0
T1 = 100.0
delta = L/32.0
lam = 0.0262
Cp = 1000
rho = 1.17 
D= lam/(rho*Cp)

print(" Heat Diffusion  ..")


#####################
def InitSolution(x):
  y=yc
  T = T0 + T1*np.exp(-0.25*((x-xc)**2+ (y- yc)**2)/delta**2)
  return T
##################
#####################
def ExactSolution(x,t):
  y=yc
  Tt = T1*np.exp(-D*t)
  T = T0 + T1*np.exp(-0.25*((x-xc)**2+ (y- yc)**2)/delta**2)
  return T
##################
def VortexLineoutX(datasets, xc, yc):
  x_out = []
  v_out = []
  T_out = []

  for d in range(len(datasets)):
    ds = yt.load(datasets[d])

    xaxis = 0  # take a line cut along x axis
    lineout = ds.ortho_ray(xaxis, (yc, 0))
    srt = np.argsort(lineout["x"])
    x = np.array(lineout["x"][srt])
    v = np.array(lineout["Ymom"][srt])
    x_out.append(x)
    v_out.append(v)
    T = np.array(lineout["temperature"][srt])
    T_out.append(T)

  return x_out, v_out, T_out
################################################

print('opening data')

# open data
list_of_files = glob.glob('./plot/*') 
latest_file = max(list_of_files, key=os.path.getctime)
print(" LAST FILE=",latest_file)
datasets = [latest_file]

#print("extracting lines at xc,yc=",xc,yc)

x, v, T  = VortexLineoutX(datasets, xc, yc)

for (xs, vs, Ts) in zip(x, v, T): 

  Te  = InitSolution(xs)
  
  # Initial  solution
  plt.plot(xs,Te,'o',label='t=0',mfc='none',markersize=3)
  plt.legend()

  # cerisse
  plt.plot(xs, Ts, label='cerisse')
  plt.ylabel("T")
  plt.xlabel(' x ', horizontalalignment='center')

  plt.show()

# Save the line plot into a file
#plot.save()


print(' .... DONE')






print(' .... DONE')
