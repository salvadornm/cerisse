import yt
import numpy as np
import matplotlib.pyplot as plt
import glob
import os

L = 1.0

xc = 0.25*L/2
yc = 0.5*L
dpdx = 0.24
rho = 1.0 
Re  = 50
h  = 0.5*L
visc = 1/Re
Vbulk = h**2*dpdx/(3*visc)

print(" Vbulk= ",Vbulk," [m/s]") 
print(" visc= ",visc," [Pa s]") 


#####################
def ExactSolution(y):
  yoh = y/L
  u = 1.5*Vbulk*(4*(1 - yoh)*yoh)  
  return u
##################
def ExtractLineoutY(datasets, xc):
  y_out = []
  v_out = []
  
  for d in range(len(datasets)):
    ds = yt.load(datasets[d])

    xaxis = 0  # take a line cut along x axis
    #lineout = ds.ortho_ray(xaxis, (yc, 0))

    yaxis = 1 # take a line cut along y axis
    lineout = ds.ortho_ray(yaxis,(xc, 0))
    

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


y, v  = ExtractLineoutY(datasets, xc)

for (ys, vs) in zip(y, v): 

  ve  = ExactSolution(ys)
  
  # Exact  solution
  plt.plot(ve,ys,'o',label='Analytical',mfc='none',markersize=3)

  # cerisse
  plt.plot(vs, ys, label='cerisse')
  plt.ylabel(' v ')
  plt.xlabel(' y/H ', horizontalalignment='center')

  plt.legend()
  plt.show()

# Save the line plot into a file
#plot.save()


print(' .... DONE')






print(' .... DONE')
