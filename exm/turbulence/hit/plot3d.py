# quick python plot to explore results from cli w/o postprocessing

import yt
import numpy as np
import matplotlib.pyplot as plt
import glob
import os
import argparse

yt.funcs.mylog.setLevel("ERROR")  # or "CRITICAL" to suppress almost everything

# take dirname as argument
parser = argparse.ArgumentParser(description="Check solution at a specific resolution.")
parser.add_argument("dir_name", type=str, help="Directory name (e.g., 'plot16')")
args = parser.parse_args()

dir_name = args.dir_name


# list of files in dir
filelist_name = f"{dir_name}/*"
list_of_files = glob.glob(filelist_name) 
latest_file = max(list_of_files, key=os.path.getctime)
print(" LAST FILE=",latest_file)

ds = yt.load(latest_file)

# Create a 2D slice plot (e.g., in the 'z' direction)
# You can change 'z' to 'x' or 'y' depending on which slice you want
slc = yt.SlicePlot(ds, 'z', ('boxlib', 'x_velocity'))

# Optional: Set contour levels or colormap
#slc.set_zlim(('boxlib', 'Density'), 0.81, 1.51)
#slc.set_cmap(('boxlib', 'Density'), 'viridis')

# Save the plot
slc.save()


