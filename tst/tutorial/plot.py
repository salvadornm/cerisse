import yt
import numpy as np
import matplotlib.pyplot as plt
#-------------------------------------------
# file to open data and plot 2D results
# and save it into a file
#-----------------------------------------

print('opening data')

# load data
ds = yt.load("plot/plt06390")

# print some information about sim
ds.print_stats()
# grid information
#g = ds.index.grids[1]
#print(g)

print('plotting ...')
print(' ##################################################### ')


# 3D
#slc = yt.SlicePlot(ds, "z", ("boxlib", "Density"))
#slc.save()

# 2D contour
p = yt.plot_2d(ds, ("boxlib", "Density"))
p.save()


print(' .... DONE')