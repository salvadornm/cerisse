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
ds = yt.load("plot/plt00050")

# list_of_files = glob.glob('./plot/*')
# latest_file = max(list_of_files, key=os.path.getctime)
# print(" LAST FILE=",latest_file)
# ds = yt.load(latest_file)

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
rho = np.array(lineout["boxlib", "temperature"][srt]) # get the variable


# plot solutions
plt.plot(np.array(lineout["index", "x"][srt]),rho,label=" cerisse ")
plt.ylabel("T")
plt.xlabel(' x ', horizontalalignment='center')


# init
# ds0 = yt.load("plot/plt00000")
# lineout = ds0.ortho_ray(xaxis, (0, 0))   # cutting through the xaxis
# srt0 = np.argsort(lineout["index", "x"]) # sort according to x position
# rho0 = np.array(lineout["boxlib", "temperature"][srt]) # get the variable
# plt.plot(np.array(lineout["index", "x"][srt]),rho0,label=" t=0 ")



plt.legend()
plt.show()

# Save the line plot into a file
#plot.save()


print(' .... DONE')
