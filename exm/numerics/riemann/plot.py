import yt
import numpy as np
import matplotlib.pyplot as plt
#-------------------------------------------
# file to open data and plot 1D results
# tst/tst1
#-----------------------------------------

print('opening data')

# load data
ds = yt.load("plot/plt00200")

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
rho = np.array(lineout["boxlib", "Density"][srt]) # get the variable

# Exact solution at t=0.2 
data = np.loadtxt('exact.dat')
plt.plot(data[:,0], data[:,1],'o',label="exact",mfc='none',markersize=3)

# plot solutions
plt.plot(np.array(lineout["index", "x"][srt]),rho,label=" cerisse ")
plt.ylabel("density")
plt.xlabel(' x ', horizontalalignment='center')


plt.legend()
plt.show()

# Save the line plot into a file
#plot.save()


print(' .... DONE')