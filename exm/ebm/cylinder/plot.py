import yt
import numpy as np
import matplotlib.pyplot as plt
import glob
import os


print(" Plotting  ..")


################################################
print('opening data')

# open data
list_of_files = glob.glob('./plot/*') 
latest_file = max(list_of_files, key=os.path.getctime)
print(" LAST FILE=",latest_file)
#datasets = [latest_file]
ds = yt.load(latest_file)
################################################
# 2D contour
p = yt.plot_2d(ds, ("boxlib", "Density"),axes_unit="unitary")
p.set_log(("boxlib", "Density"), False)
p.save()


################################################

print(' .... DONE')
