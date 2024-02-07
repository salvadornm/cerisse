import yt
import numpy as np
import matplotlib.pyplot as plt

##################
def VortexLineout(datasets, xc, yc):
  x_out = []
  v_out = []
  y_out = []
  u_out = []

  for d in range(len(datasets)):
    ds = yt.load(datasets[d])

    xaxis = 0  # take a line cut along x axis
    lineout = ds.ortho_ray(xaxis, (yc, 0))
    srt = np.argsort(lineout["x"])
    x = np.array(lineout["x"][srt])
    v = np.array(lineout["velocity_y"][srt])
    x_out.append(x)
    v_out.append(v)

    xaxis = 1  # take a line cut along y axis
    lineout = ds.ortho_ray(xaxis, (0, xc))
    srt = np.argsort(lineout["y"])
    y = np.array(lineout["y"][srt])
    u = np.array(lineout["velocity_x"][srt])
    y_out.append(y)
    u_out.append(u)

  return x_out, v_out, y_out, u_out
################################################

v0 = 0.05 * np.sqrt(1.4 * 287 * 300) * 100 # convection velocity 1735.95
beta = 0.02  # vortex strength
R = 0.5
#xc = (5.0 + ((0.288 * v0) % 10.0)) % 10.0
xc = 5.0
yc = 5.0

print('opening data')

# open data
datasets = ["plt04500"]

print('extracting lines at xc,yc=',xc,yc)


x, v, y, u = VortexLineout(datasets, xc, yc)

xp = np.linspace(0, 10, 300)
u_exact = -v0*beta*(xp-yc)/R*np.exp(-0.5*((xp-yc)/R)**2)
v_exact = v0*beta*(xp-xc)/R*np.exp(-0.5*((xp-xc)/R)**2)

# print('x=',x)
# print('v=',v)

for (xs, vs) in zip(x, v):
  plt.plot(xs, vs, label='AMR 64*2')

# exact solution
plt.plot(xp,v_exact,'o',label='exact')
plt.legend()
plt.show()






print(' .... DONE')
