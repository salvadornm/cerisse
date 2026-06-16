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
ds_teno = yt.load("plotTENO5")
ds_skew = yt.load("plotSKEW4")
ds_weno5fs = yt.load("plotWENOZ5")

ds_godunov = yt.load("plotN_GODUNOV")
ds_plm = yt.load("plotN_PLM")
ds_weno3 = yt.load("plotN_WENO3")
ds_weno5 = yt.load("plotN_WENOJS5")
ds_wenoz5 = yt.load("plotN_WENOZ5")
ds_teno5 = yt.load("plotN_TENO5")




# print some information about sim
ds.print_stats()
# grid information
g = ds.index.grids[1]

print(g)

print('plotting ...')
print(' ##################################################### ')


xaxis = 0  # take a line cut along the x axis
lineout = ds.ortho_ray(xaxis, (0, 0))   # cutting through the xaxis

lineout_teno = ds_teno.ortho_ray(xaxis, (0, 0))   # cutting through the xaxis
lineout_skew = ds_skew.ortho_ray(xaxis, (0, 0))   # cutting through the xaxis       
lineout_godunov = ds_godunov.ortho_ray(xaxis, (0, 0))   # cutting through the xaxis       
lineout_plm = ds_plm.ortho_ray(xaxis, (0, 0))   # cutting through the xaxis       
lineout_weno3 = ds_weno3.ortho_ray(xaxis, (0, 0))   # cutting through the xaxis       
lineout_weno5 = ds_weno5.ortho_ray(xaxis, (0, 0))   # cutting through the xaxis       
lineout_wenoz5 = ds_wenoz5.ortho_ray(xaxis, (0, 0))   # cutting through the xaxis       
lineout_teno5 = ds_teno5.ortho_ray(xaxis, (0, 0))   # cutting through the xaxis   
lineout_weno5fs = ds_weno5fs.ortho_ray(xaxis, (0, 0))   # cutting through the xaxis

srt = np.argsort(lineout["index", "x"]) # sort according to x position
srt_teno = np.argsort(lineout_teno["index", "x"]) # sort according to x position
srt_skew = np.argsort(lineout_skew["index", "x"]) # sort according to x position
srt_godunov = np.argsort(lineout_godunov["index", "x"]) # sort according to x position
srt_plm = np.argsort(lineout_plm["index", "x"]) # sort according to x position
srt_weno3 = np.argsort(lineout_weno3["index", "x"]) # sort according to x position
srt_weno5 = np.argsort(lineout_weno5["index", "x"]) # sort according to x position
srt_wenoz5 = np.argsort(lineout_wenoz5["index", "x"]) # sort according to x position
srt_teno5 = np.argsort(lineout_teno5["index", "x"]) # sort according to x position
srt_weno5fs = np.argsort(lineout_weno5fs["index", "x"]) # sort according to x position

rho = np.array(lineout["boxlib", "Density"][srt]) # get the variable
rho_teno = np.array(lineout_teno["boxlib", "Density"][srt_teno]) # get the variable
rho_skew = np.array(lineout_skew["boxlib", "Density"][srt_skew]) # get the variable
rho_godunov = np.array(lineout_godunov["boxlib", "Density"][srt_godunov]) # get the variable
rho_plm = np.array(lineout_plm["boxlib", "Density"][srt_plm]) # get the variable
rho_weno3 = np.array(lineout_weno3["boxlib", "Density"][srt_weno3]) # get the variable
rho_weno5 = np.array(lineout_weno5["boxlib", "Density"][srt_weno5]) # get the variable
rho_wenoz5 = np.array(lineout_wenoz5["boxlib", "Density"][srt_wenoz5]) # get the variable
rho_teno5 = np.array(lineout_teno5["boxlib", "Density"][srt_teno5]) # get the variable
rho_weno5fs = np.array(lineout_weno5fs["boxlib", "Density"][srt_weno5fs]) # get the variable

# Exact solution at t=0.2 
data = np.loadtxt('exact.dat')
plt.plot(data[:,0], data[:,1],'o',label="exact",mfc='none',markersize=4)

# plot solutions
#plt.plot(np.array(lineout["index", "x"][srt]),rho,label=" cerisse ")
#plt.plot(np.array(lineout_teno["index", "x"][srt_teno]),rho_teno,label=" TENO ")
#plt.plot(np.array(lineout_skew["index", "x"][srt_skew]),rho_skew,label=" SKEW ")
#plt.plot(np.array(lineout_godunov["index", "x"][srt_godunov]),rho_godunov,label=" Godunov ")
#plt.plot(np.array(lineout_plm["index", "x"][srt_plm]),rho_plm,label=" PLM ")
#plt.plot(np.array(lineout_weno3["index", "x"][srt_weno3]),rho_weno3,label=" WENO3 ")
#plt.plot(np.array(lineout_weno5["index", "x"][srt_weno5]),rho_weno5,label=" WENOJS5 ")
plt.plot(np.array(lineout_wenoz5["index", "x"][srt_wenoz5]),rho_wenoz5,label=" Recons WENOZ5 ")
#plt.plot(np.array(lineout_teno5["index", "x"][srt_teno5]),rho_teno5,label=" TENO5 ")
plt.plot(np.array(lineout_weno5fs["index", "x"][srt_weno5fs]),rho_weno5fs,label=" Flux Split WENOZ5 ")

plt.ylabel("density")
plt.xlabel(' x ', horizontalalignment='center')


plt.legend()
plt.show()

# Save the line plot into a file
#plot.save()


print(' .... DONE')