# Tutorial

This  page will explain in detail the set-up, run and visualization of a typical case

As a glace, domain dimension and control parameters are handled in file **inputs**
while proiblem description is in **prob.h**


## Set-up

The problem involves a heavy fluid falling into a light fluid.
The upper half of the domain is filled with a fluid of density 2, while the lower part is filled witha fluid of density 1. The initial pressure distrbution follows hydrostatic and a velocity perturbation initiates the instability.

TABLE INITIAL CONDITIONS

## prob.h


### Set-up of Problem Parameters

The parameters of the problem (not the domain) are wrapped into a **ProbParm** structure

```cpp
// problem parameters
struct ProbParm {
  Real p_1 = 2.0;
  Real p_2 = 2.0;
  Real rho_1 = 1.0;
  Real rho_2 = 2.0;
  Real tau = 6.0;
  Real M0 = 0.1;
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
  Real Pint = prob_parm.p_2 -prob_parm.rho_2*prob_parm.grav*Lint; // interface Pressure
  // bottom
  if (y < Lint) {
    rhot = prob_parm.rho_1;
    Pt = Pint - rhot*prob_parm.grav*y;
    uxt = 0;
  } 
  // top
  else {
    rhot = prob_parm.rho_2;
    Pt = prob_parm.p_2 - rhot*prob_parm.grav*y;
    uxt = 0;
  }
  state(i, j, k, cls.URHO) = rhot;
  state(i, j, k, cls.UMX) = rhot * uxt;
  state(i, j, k, cls.UMY) = rhot * uyt;
  Real et = Pt / (cls.gamma - Real(1.0));
  state(i, j, k, cls.UET) = et + Real(0.5) * rhot * (uxt * uxt + uyt * uyt); 
```


The array state holds

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

### Paraview





