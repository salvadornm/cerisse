import yt
import numpy as np
import matplotlib.pyplot as plt
import glob
import os


print(" Plotting  EB ..")


################################################
print('opening data')

# open data
list_of_files = glob.glob('./plot/*') 
latest_file = max(list_of_files, key=os.path.getctime)
print(" LAST FILE=",latest_file)
ds = yt.load(latest_file)

# Define a new derived field that masks Density where vfrac < 0.5
def masked_density(field, data):
    vfrac = data[("boxlib", "vfrac")]
    dens = data[("boxlib", "Density")]
    return np.where(vfrac >= 0.1, dens, np.nan)


ds.add_field(
    ("boxlib", "masked_density"),
    function=masked_density,
    sampling_type="cell",
    force_override=True,
)

################################################
# 2D contour (will remove solid points)
p = yt.plot_2d(ds, ("boxlib", "masked_density"),axes_unit="unitary")
p.set_log(("boxlib", "masked_density"), False)
result = p.save()

# image_path = result[0][1]  # This gets the full path to the saved image
# # Load the image and show it using matplotlib
# import matplotlib.image as mpimg
# img = mpimg.imread(image_path)

# plt.imshow(img)
# plt.axis('off')
# plt.title("Masked Density Plot")
# plt.show()


################################################

print(' .... DONE')
