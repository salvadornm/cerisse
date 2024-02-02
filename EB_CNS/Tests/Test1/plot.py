import yt
import numpy as np
import matplotlib.pyplot as plt
#-------------------------------------------
# file to open data and plot 1D results
# Examples/Test1
#-----------------------------------------

print('opening data')

# load data
ds = yt.load("plt00200")

# print some information about sim
ds.print_stats()
# grid information
g = ds.index.grids[1]

print(g)

print('plotting ...')
print(' ##################################################### ')


xaxis = 0  # take a line cut along the x axis
lineout = ds.ortho_ray(xaxis, (0, 0))   # cutting through the xaxis


srt = np.argsort(lineout["index", "x"]) # sort according to x position
rho = np.array(lineout["gas", "density"][srt]) # get the variable

# plot solutions
plt.plot(np.array(lineout["index", "x"][srt]),rho,label=" WENO ")
plt.ylabel("density")
plt.xlabel(' x ', horizontalalignment='center')

# Exact solution at t=0.2  columns c(0)=x c(1)=dens c(2)=pres c(3)=velx c(4)=eint c(5)-csd
data = np.loadtxt('exact.dat')
plt.plot(data[:,0], data[:,1],'ro',label="exact")

plt.legend()
plt.show()

# Save the line plot into a file
#plot.save()


print(' .... DONE')