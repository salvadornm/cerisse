#!/usr/bin/env python

# Import modules
import argparse
import sys
import numpy as np
import scipy.interpolate as spi
import matplotlib.pyplot as plt

# Some defaults variables
# plt.rc("text", usetex=True)
# plt.rc("font", family="serif", serif="Times")
plt.rc("font", size=14)
cmap_med = ["#F15A60", "#7AC36A", "#5A9BD4", "#FAA75B", "#9E67AB", "#CE7058", "#D77FB4", "#737373"]
cmap = ["#EE2E2F", "#008C48", "#185AA9", "#F47D23", "#662C91", "#A21D21", "#B43894", "#010202"]
dashseq = [(None), [10, 5], [10, 4, 3, 4], [3, 3], [10, 4, 3, 4, 3, 4], [3, 3], [3, 3]]
markertype = ["s", "d", "o", "p", "h"]

# Some helper functions
def div0(a, b):
  """ Ignore division by 0, just replace it by 0,
  From: http://stackoverflow.com/questions/26248654/numpy-return-0-with-divide-by-zero
  e.g. div0( [-1, 0, 1], 0 ) -> [0, 0, 0]
  """
  with np.errstate(divide="ignore", invalid="ignore"):
    c = np.true_divide(a, b)
    c[~np.isfinite(c)] = 0  # -inf inf NaN
  return c

def abs2(x):
  """This is equivalent to np.abs(x)**2 or x*np.conj(x)
  To make it faster, add this right before the function definition
  import numba
  @numba.vectorize([numba.float64(numba.complex128),numba.float32(numba.complex64)])
  """
  return x.real**2 + x.imag**2

def vel2Ek(uf, vf, wf):
  """Get energy spectrum from velocity field. uf, vf, wf are Fourier transformed velocities:
  uf = np.fft.rfftn(u, s=(N, N, N))
  vf = np.fft.rfftn(v, s=(N, N, N))
  wf = np.fft.rfftn(w, s=(N, N, N))
  """ 
  N = uf.shape[0]

  k = np.concatenate((np.arange(int(round(0.5 * N))), np.arange(-int(round(0.5 * N)), 0, 1)), axis=0)
  khalf = np.arange(int(round(0.5 * N)) + 1)
  k1, k2, k3 = np.meshgrid(k, k, khalf, indexing="ij")
  kmag = np.sqrt(k1 ** 2 + k2 ** 2 + k3 ** 2)

  kbins = np.arange(1, int(round(0.5 * N)) + 1)
  Nbins = len(kbins)
  whichbin = np.digitize(kmag.flat, kbins)
  ncount = np.bincount(whichbin)

  KI = (abs2(uf) + abs2(vf) + abs2(wf)) * 0.5 / N ** 6
  KI[:, :, 1:-1] += (abs2(uf[:, :, 1:-1]) + abs2(vf[:, :, 1:-1]) + abs2(wf[:, :, 1:-1])) * 0.5 / N ** 6

  Eku = np.zeros(len(ncount) - 1)
  for n in range(1, len(ncount)):
    Eku[n - 1] = np.sum(KI.flat[whichbin == n])

  ku = 0.5 * (kbins[0 : Nbins - 1] + kbins[1:Nbins]) + 1
  Eku = Eku[1:Nbins]
  return ku, Eku

################################### Start of the program ###################################
# Parse arguments
parser = argparse.ArgumentParser(description="Generate the velocity fluctuations for the HIT IC")
parser.add_argument("-k0", help="Wave number containing highest energy", type=float, default=4.0)
parser.add_argument("-N", help="Output resolution", type=int, default=32)
parser.add_argument("-m", help="Generation resolution / output resolution", type=int, default=8)
parser.add_argument("-s", "--seed", help="Random number generator seed", type=int, default=42)
parser.add_argument("-p", "--plot", help="Save plots", action="store_true")
args = parser.parse_args()

xlo = 0
xhi = 2.0 * np.pi

# Define spectrum function
# spectrum = lambda k: (
#       16.0 * np.sqrt(2.0 / np.pi)
#       * (k ** 4) / (args.k0 ** 5)
#       * np.exp(-2.0 * (k ** 2) / (args.k0 ** 2)))

eta = 45 # Kolmogorov wavenumber
spectrum = lambda k: (3.07 / np.sqrt(np.pi) * (12/5)**2 
                      * (k / args.k0)**4 / (1 + (k / args.k0 * np.sqrt(12/5))**2)**(17/6) 
                      * np.exp(-2.0 * (k/eta)**2))

# 1. Generate random velocity on a big cube 
N = args.N * args.m
halfN = int(round(0.5 * N))
L = xhi - xlo
dx = L / N

# Only work if N and args.N are even
if not ((args.N % 2 == 0) and N % 2 == 0):
  print("N or args.N is not even. Exiting")
  sys.exit(1)

# Get cell centered values and meshed grid
x = np.linspace(xlo, xhi, N + 1)
xc = (x[1:] + x[:-1]) / 2  # get cell center coordinates
X, Y, Z = np.meshgrid(xc, xc, xc, indexing="ij")

# Get the wave numbers and associated quantities
k = np.concatenate((np.arange(halfN), np.arange(-halfN, 0, 1)), axis=0)
khalf = np.arange(halfN + 1)
k1, k2, k3 = np.meshgrid(k, k, khalf, indexing="ij")
kmag = np.sqrt(k1 ** 2 + k2 ** 2 + k3 ** 2)
k12 = np.sqrt(k1 ** 2 + k2 ** 2)
k1k12 = div0(k1, k12)
k2k12 = div0(k2, k12)
k3kmag = div0(k3, kmag)
k12kmag = div0(k12, kmag)

# Generate data
# Energy spectrum
Ek = spectrum(kmag)

# Draw random numbers
np.random.seed(args.seed)
phi1 = np.random.uniform(0, 2 * np.pi, np.shape(kmag))
phi2 = np.random.uniform(0, 2 * np.pi, np.shape(kmag))
phi3 = np.random.uniform(0, 2 * np.pi, np.shape(kmag))

# the random quantities
prefix = np.sqrt(2.0 * div0(Ek, 4.0 * np.pi * (kmag ** 2)))
a = prefix * np.exp(1j * phi1) * np.cos(phi3)
b = prefix * np.exp(1j * phi2) * np.sin(phi3)

# the random velocities
uf = k2k12 * a + k1k12 * k3kmag * b
vf = k2k12 * k3kmag * b - k1k12 * a
wf = -k12kmag * b

# Impose the 3D spherical symmetry (to ensure we have a real signal)
# equiv: uf[-l,-m,0] = np.conj(uf[ l, m,0]) for l=0..N/2 and m=0..N/2
uf[N:halfN:-1, N:halfN:-1, 0] = np.conj(uf[1:halfN, 1:halfN, 0])
# symmetry on first column
uf[N:halfN:-1, 0, 0] = np.conj(uf[1:halfN, 0, 0])
# symmetry on first row
uf[0, N:halfN:-1, 0] = np.conj(uf[0, 1:halfN, 0])
# symmetry about the (halfN,halfN) element
uf[halfN - 1 : 0 : -1, N : halfN - 1 : -1, 0] = np.conj(
    uf[halfN + 1 : N, 1 : halfN + 1, 0])

vf[N:halfN:-1, N:halfN:-1, 0] = np.conj(vf[1:halfN, 1:halfN, 0])
vf[halfN - 1 : 0 : -1, N : halfN - 1 : -1, 0] = np.conj(
    vf[halfN + 1 : N, 1 : halfN + 1, 0])
vf[N:halfN:-1, 0, 0] = np.conj(vf[1:halfN, 0, 0])
vf[0, N:halfN:-1, 0] = np.conj(vf[0, 1:halfN, 0])

wf[N:halfN:-1, N:halfN:-1, 0] = np.conj(wf[1:halfN, 1:halfN, 0])
wf[halfN - 1 : 0 : -1, N : halfN - 1 : -1, 0] = np.conj(
    wf[halfN + 1 : N, 1 : halfN + 1, 0])
wf[N:halfN:-1, 0, 0] = np.conj(wf[1:halfN, 0, 0])
wf[0, N:halfN:-1, 0] = np.conj(wf[0, 1:halfN, 0])

# Normalize. Because we are generating the data in wavenumber space,
# we have to multiply by N**3 because in the definition of the numpy
# ifftn there is a 1/N**n.
uf = uf * N ** 3
vf = vf * N ** 3
wf = wf * N ** 3

# Quick check on energy content (make sure you add both the current
# contribution and the one we are neglecting because we are assuming
# real input data)
print('Energy = ∫ E(k) dk = 0.5 * ∫ (uf**2 + vf**2 wf**2) dk1 dk2 dk3 = {0:.10f} ~= 3/2'.format(
    (np.sum(abs2(uf          ) + abs2(vf          ) + abs2(wf          )) +
     np.sum(abs2(uf[:,:,1:-1]) + abs2(vf[:,:,1:-1]) + abs2(wf[:,:,1:-1])))
    * 0.5 / N**6))

# 2. Filter to the desired resolution
uo = np.fft.irfftn(uf, s=(N, N, N))
vo = np.fft.irfftn(vf, s=(N, N, N))
wo = np.fft.irfftn(wf, s=(N, N, N))

# box filter
fac = args.m
ua = np.zeros((args.N, args.N, args.N))
va = np.zeros((args.N, args.N, args.N))
wa = np.zeros((args.N, args.N, args.N))
ksgs = np.zeros_like(ua)
for i in range(args.N):
  for j in range(args.N):
    for k in range(args.N):
      ua[i, j, k] = np.mean(uo[fac*i:fac*(i+1), fac*j:fac*(j+1), fac*k:fac*(k+1)])
      va[i, j, k] = np.mean(vo[fac*i:fac*(i+1), fac*j:fac*(j+1), fac*k:fac*(k+1)])
      wa[i, j, k] = np.mean(wo[fac*i:fac*(i+1), fac*j:fac*(j+1), fac*k:fac*(k+1)])
      uaf = uo[fac*i:fac*(i+1), fac*j:fac*(j+1), fac*k:fac*(k+1)] - ua[i, j, k]
      vaf = vo[fac*i:fac*(i+1), fac*j:fac*(j+1), fac*k:fac*(k+1)] - va[i, j, k]
      waf = wo[fac*i:fac*(i+1), fac*j:fac*(j+1), fac*k:fac*(k+1)] - wa[i, j, k]
      ksgs[i, j, k] = 0.5 * np.mean(uaf**2 + vaf**2 + waf**2)

# check energy
kea = 0.5 * (ua**2 + va**2 + wa**2)
print("Energy res   = {0:.10f}".format(np.sum(kea)*(1 / args.N)**3))
print("Energy sgs   = {0:.10f}".format(np.sum(ksgs)*(1 / args.N)**3))
print("Energy total = {0:.10f} ~= 3/2".format(np.sum(kea + ksgs)*(1 / args.N)**3))
# if we use interp downsample u, the total ke will be > 3/2

xr = np.linspace(xlo, xhi, args.N + 1)
xrc = (xr[1:] + xr[:-1]) / 2
Xr, Yr, Zr = np.meshgrid(xrc, xrc, xrc, indexing="ij")

# Save the data in Fortran ordering
fname = "hit2_ic_{0:d}_{1:d}.dat".format(int(args.k0), args.N)
data = np.vstack((Xr.reshape(-1, order="F"), Yr.reshape(-1, order="F"), Zr.reshape(-1, order="F"),
                  ua.reshape(-1, order="F"), va.reshape(-1, order="F"), wa.reshape(-1, order="F"),
                  ksgs.reshape(-1, order="F"))).T
np.savetxt(fname, data, fmt="%.18e", delimiter=",", header="x, y, z, u, v, w, ksgs")

if args.plot:
  # plot velocity field
  datmin = uo.min()
  datmax = uo.max()
  def myimshow(ax, u, title):
    ax.imshow(u, origin="lower", extent=[xlo, xhi, xlo, xhi],
              cmap="RdBu_r", vmin=datmin, vmax=datmax)
    ax.set_title(title)

  fig, ax = plt.subplots(nrows=3, ncols=3, figsize=(9, 10), layout="constrained")
  myimshow(ax[0, 0], uo[:, :, 0].T, r"Original $u (x,y)$")
  myimshow(ax[1, 0], uo[:, 0, :].T, r"Original $u (x,z)$")
  myimshow(ax[2, 0], uo[0, :, :].T, r"Original $u (y,z)$")
  myimshow(ax[0, 1], ua[:, :, 0].T, r"Filtered $u (x,y)$")
  myimshow(ax[1, 1], ua[:, 0, :].T, r"Filtered $u (x,z)$")
  myimshow(ax[2, 1], ua[0, :, :].T, r"Filtered $u (y,z)$")
  myimshow(ax[0, 2], np.sqrt(2*ksgs[:, :, 0].T), r"$\sqrt{2k_{sgs}} (x,y)$")
  myimshow(ax[1, 2], np.sqrt(2*ksgs[:, 0, :].T), r"$\sqrt{2k_{sgs}} (x,z)$")
  myimshow(ax[2, 2], np.sqrt(2*ksgs[0, :, :].T), r"$\sqrt{2k_{sgs}} (y,z)$")

  img_name = "hit2_ic_u_{0:d}_{1:d}.png".format(int(args.k0), args.N)
  print("Writing to " + img_name + "...")
  plt.savefig(img_name, format="png")

  # plot spectra
  Eko = spectrum(khalf)
  ku, Eku = vel2Ek(uf, vf, wf)
  ui2 = np.fft.rfftn(ua, s=(args.N, args.N, args.N))
  vi2 = np.fft.rfftn(va, s=(args.N, args.N, args.N))
  wi2 = np.fft.rfftn(wa, s=(args.N, args.N, args.N))
  ka, Eka = vel2Ek(ui2, vi2, wi2)

  plt.figure()
  plt.loglog(khalf, Eko, linestyle=(0,dashseq[0]), color=cmap[-1], lw=2, label="Analytical")
  plt.loglog(ku, Eku, linestyle=(0,dashseq[1]), color=cmap[0], lw=2, label="Unfiltered")
  plt.loglog(ka, Eka, linestyle=(0,dashseq[2]), color=cmap[1], lw=2, label="Filtered")
  # plt.xlim([None, 50])
  plt.ylim([1e-8, 10])
  plt.xlabel(r"$k$", fontsize=22, fontweight="bold")
  plt.ylabel(r"$E(k)$", fontsize=22, fontweight="bold")
  # plt.setp(ax.get_xmajorticklabels(), fontsize=18, fontweight="bold")
  # plt.setp(ax.get_ymajorticklabels(), fontsize=18, fontweight="bold")
  plt.legend()

  img_name = "hit2_ic_spectrum_{0:d}_{1:d}.png".format(int(args.k0), args.N)
  print("Writing to " + img_name + "...")
  plt.savefig(img_name, format="png")

print("Done!")