# Tutorial

This  page will explain in detail the set-up, run and visualization of a typical case

As a glace, domain dimension and control parameters are handled in file **inputs**
while problem description is in **prob.h**


## Set-up

The problem involves a heavy fluid falling into a light fluid.
The upper half of the domain is filled with a fluid of density 2, while the lower part is filled with a fluid of density 1. The initial pressure distribution follows hydrostatic and a velocity perturbation initiates the instability.

TABLE INITIAL CONDITIONS

CFL=0.3, time = 2

Initial conditions. Following Shi et al [^1]


Meshes 64x256 and 128x512 

<img src="../../images/tutorial_RT1.png" width=300 height=700>
<img src="../../images/tutorial_RT2.png" width=300 height=700>


## prob.h


### Set-up of Problem Parameters

The parameters of the problem (not the domain) are wrapped into a **ProbParm** structure

```cpp
// problem parameters
struct ProbParm {
  Real p_int = 2.0;
  Real rho_1 = 1.0;
  Real rho_2 = 2.0;
  Real grav = -1.0; 
  Real eps =  0.025;
};
```

### Set-up of  initial conditions

This is done in the function **prob_initdata**

```cpp
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
prob_initdata(int i, int j, int k, Array4<Real> const &state,
              GeometryData const &geomdata, ProbClosures const &cls,
              ProbParm const &prob_parm) {
  const Real *prob_lo = geomdata.ProbLo();
  const Real *prob_hi = geomdata.ProbHi();
  const Real *dx = geomdata.CellSize();

  Real x = prob_lo[0] + (i + Real(0.5)) * dx[0];
  Real y = prob_lo[1] + (j + Real(0.5)) * dx[1];
```
where the spatial coordinates **x** and **y** of the cells are determined from the mesh sizes and problem dimensions (which are defined in **inputs**).

The initial conditions are defined as

```cpp
  Real Pt, rhot, uxt,uyt;
  Real Lint = prob_hi[0] / 2;
  Real Pint = prob_parm.p_int; // interface Pressure

  const Real freq = Real(8)*Real(3.14159265359); // wavelength = x-domain

  Real yrel = y - Lint;
  Real delta= 0.2*Lint;   // region size where perturbation is significant
  Real delta2  = dx[1]/5; // transition region between top/bottom
  Real step = Real(0.5) + Real(0.5)*tanh(yrel/delta2);
  rhot = step*prob_parm.rho_2 + (Real(1.0) -step)*prob_parm.rho_1;
  Pt = Pint + rhot*prob_parm.grav*(y - Lint); // hydrostatic pressure

  uxt = Real(0.0);
  uyt = -prob_parm.eps*cos(freq*x)*aux; // perturbation in y-component

  state(i, j, k, cls.URHO) = rhot;
  state(i, j, k, cls.UMX) = rhot * uxt;
  state(i, j, k, cls.UMY) = rhot * uyt;
  Real et = Pt / (cls.gamma - Real(1.0));
  state(i, j, k, cls.UET) = et + Real(0.5) * rhot * (uxt * uxt + uyt * uyt); 
```


The array **state** holds all the information required

### Solvers


### Source


### AMR


## Compile and Running


## Post-processing

### Python



### Visit

In this example, we will use Visit (recall that Visit cannot be used in 1D).
To load the data, open header files of a particular directory

FIGURE HERE


If the instant 0 is open, the results for density will look like 

FIGURE HERE


Alternatively you can use the script **cerisse**, to open all directories at the same time (to make an animation for example). See  [Tips](tips.md) to set-up the script.

``` bash
$ cerisse visit
```

<img src="../../images/tutorial_RT_amr2.png" width=400 height=400>
<img src="../../images/tutorial_RT_amr3.png" width=400 height=400>



### Paraview



[^1]:  J Shi et al, *Resolution of high order WENO schemes for complicated flow structures*, J Computational Physics, (2003), 186, pp 690-696.
https://www.sciencedirect.com/science/article/pii/S0021999103000949
