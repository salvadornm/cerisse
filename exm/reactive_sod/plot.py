import numpy as np
import matplotlib.pyplot as plt
import yt
import glob

folders = ["muscl", "weno", "teno"]

# Plot the results
time_int = [16, 23]
var_name = ["temperature", "pressure", "rho_AR", "rho_OH"]
line_styles = ["r-", "b--", "g-."]

plt.figure(figsize=(6, 8))
plt.suptitle("$t = " + str([t * 10 for t in time_int]) + "\, \mu s$")
for j, folder in enumerate(folders):
    fnames = sorted(glob.glob(folder + "/plt*"))
    for t in time_int:
        ds = yt.load(fnames[t])
        xaxis = 0  # take a line cut along the x axis
        lineout = ds.ortho_ray(xaxis, (0, 0))  # cutting through the xaxis
        srt = np.argsort(lineout["index", "x"])  # sort according to x position

        # plot solutions
        for i, name in enumerate(var_name):
            plt.subplot(len(var_name), 1, i + 1)
            (line,) = plt.plot(
                np.array(lineout["index", "x"][srt]),
                np.array(lineout["boxlib", name][srt]),
                line_styles[j],
            )
            plt.ylabel(name)
            plt.xlim(
                [float(ds.domain_left_edge[0]), float(ds.domain_right_edge[0])]
            )
    line.set_label(folder)

# Plot reference data
# note: a slightly different mechanism is used in ref
ref_data = np.loadtxt("ref_data", skiprows=1)
plt.subplot(len(var_name), 1, 1 + var_name.index("temperature"))
plt.plot(ref_data[:, 0], ref_data[:, 1], "ko")
plt.plot(ref_data[:, 0], ref_data[:, 4], "ko")
plt.subplot(len(var_name), 1, 1 + var_name.index("pressure"))
plt.plot(ref_data[:, 0], ref_data[:, 2], "ko")
plt.plot(ref_data[:, 0], ref_data[:, 5], "ko")
plt.subplot(len(var_name), 1, 1 + var_name.index("rho_AR"))
plt.plot(ref_data[:, 0], ref_data[:, 3] * 0.88586543, "ko")
plt.plot(ref_data[:, 0], ref_data[:, 6] * 0.88586543, "ko")

plt.subplot(len(var_name), 1, len(var_name))
plt.plot([-1], [0], "o", color="k", label="ref")
plt.legend()
plt.xlabel("x [m]")
plt.tight_layout()

plt.savefig("compare_ref.png", facecolor="white")
print("Saved plot to compare_ref.png!")
