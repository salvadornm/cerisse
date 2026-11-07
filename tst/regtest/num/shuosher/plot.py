import yt
import numpy as np
import matplotlib.pyplot as plt
import glob
import os

#-------------------------------------------
# file to open data and plot 1D results
#-----------------------------------------

print('opening data')

# load data
#ds = yt.load("plot/plt00200")

# open data
list_of_files = glob.glob('./plot/*')
latest_file = max(list_of_files, key=os.path.getctime)
print(" LAST FILE=",latest_file)
datasets = [latest_file]
ds = yt.load(datasets[0])

# print some information about sim
ds.print_stats()
# grid information
g = ds.index.grids[1]

print(g)

print('plotting ...')
print(' ##################################################### ')


# load reference data plotref
ds0 = yt.load("plotREF")

xaxis = 0  # take a line cut along the x axis
lineout = ds.ortho_ray(xaxis, (0, 0))   # cutting through the xaxis
lineout0 = ds0.ortho_ray(xaxis, (0, 0))   # cutting through the xaxis REF


srt = np.argsort(lineout["index", "x"]) # sort according to x position
rho = np.array(lineout["boxlib", "Density"][srt]) # get the variable

# Reference solution at t=1.8 
srt0 = np.argsort(lineout0["index", "x"]) # sort according to x position
rho0 = np.array(lineout0["boxlib", "Density"][srt0]) # get the variable



# plot solutions
plt.plot(np.array(lineout0["index", "x"][srt0]),rho0,label=" REF ")


plt.plot(np.array(lineout["index", "x"][srt]),rho,label=" cerisse ")

plt.ylabel("density")
plt.xlabel(' x ', horizontalalignment='center')


plt.legend()
plt.show()

# Save the line plot into a file
#plot.save()


print(' .... DONE')
