import yt
import numpy as np
import matplotlib.pyplot as plt
#-------------------------------------------
# file to open data and plot 1D results
# Examples/Test1
#-----------------------------------------

print('opening data')

# load data
ds = yt.load("plt00100")

# print some information about sim
ds.print_stats()
# grid infi
g = ds.index.grids[1]
print(g)

print('plotting ...')


#plot = yt.LinePlot(
#  ds, [("gas", "density")],(0.0, 0.0, 0.0), (1.0, 0.0, 0.0),400 )

#profiles = []
# Create a data container to hold the whole dataset.
#ad = ds.all_data()
#profiles.append(
#  yt.create_profile(ad, [("gas", "x")],fields=[("gas", "density")])
#)
#plt = yt.ProfilePlot.from_profiles(profiles)


#plt.xlabel(' x ', horizontalalignment='center')

# Save the line plot into a file
#plot.save()




print(' .... DONE')