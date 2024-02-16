# Running


NEW CODE DISCUSSION DOCUMENT


Folde


## Building the code

There are two ways to compile the code, using **GNUMake** or **CMake**. The former is recommended because it is simple and more comprehensible, while the latter is mainly used for automated testing with CTest.

## Required files

Each case must contain the following 5 files to get compiled and run

* **`GNUmakefile`** - sets the compile-time options
* **`prob.H`** - define functions for initialization (`prob_initdata`), boundary conditions (`bcnormal`), etc.
* **`prob.cpp`** - define function to initialize AMReX (`amrex_probinit`)
* **`prob_parm.H`** - define the `ProbParm` struct, which contains data for initialization or boundary conditions
* **`inputs`** - contain runtime options for the cerisse executable. Some may contain the `.ini` extension, but they are the same thing

There may also be a few optional files

* **`Make.package`** - link the `.H` and `.cpp` files to the compiler. It can be absorbed into `GNUmakefile`
* **`CMakeLists.txt`** - for building with CMake

We will now go through these files one by one.

## GNUmakefile

This file sets general options, usually modified once and that are required to compile the code. It is divided into sections

**AMReX** options

```makefile
# AMReX
DIM = 1             # dimension of the problem (1/2/3)
COMP = gnu          # compiler (gnu/intel/..)
PRECISION = DOUBLE  # floating-point precision (it has not been tested in SINGLE)
```

**PATH**


These variable select the path and had to be local

```makefile
AMR_SOLVER ?= $(abspath ../../)
AMREX_HOME ?= $(abspath ../../lib/amrex)
```

## Problem Set-up example


```
#ifndef CNS_PROB_H_
#define CNS_PROB_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_ParmParse.H>
#include <Closures.h>
#include <RHS.h>

using namespace amrex;
```
Start of the problem, these line are always the same


```
/ problem parameters
struct ProbParm {
  Real p_l = 1.0;
  Real p_r = 0.1;
  Real rho_l = 1.0;
  Real rho_r = 0.125;
  Real u_l = 0.0;
  Real u_r = 0.0;
};
```

Definition of problem parameters that may be useful



```
inline Vector<std::string> cons_vars_names={"Density","Xmom","Ymom","Zmom","Energy"};
```

Vector of Name of variables to solve



```
inline Vector<int> cons_vars_type={0,1,2,3,0};
```

Not known ??



```
typedef closures_dt<indicies_t, visc_suth_t, cond_suth_t,
                    calorifically_perfect_gas_t<indicies_t>>
    ProbClosures;
typedef rhs_dt<riemann_t<false, ProbClosures>, no_diffusive_t, no_source_t>
    ProbRHS;
```

These line define the problem, by firts defining the closures

### Closures


| Closure                      | Type          | Use | Description                                                  |
| --------------------------- | ------------- |:-------:| ------------------------------------------------------------ |
| ```visc_suth```             |            |   often      | Sutherland Viscosity                           |
| ```cond_suth```             |            |   often      | Sutherland Conductivty                         |
| ```calorifically_perfect_gas_t```             |            |   often      | perfect gas                      |


### RHS

| RHS                     | Type          | Use | Description                                                  |
| --------------------------- | ------------- |:-------:| ------------------------------------------------------------ |
| ```riemann_t```             |            |   often      | Rimeann Problem                          |
| ```no_diffusive_t```        |            |   often      | No diffusive part (Euler)                |
| ```no_source_t```           |            |   often      | No source term    |



## Problem definition

```
// initial condition
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
prob_initdata(int i, int j, int k, Array4<Real> const &state,
              GeometryData const &geomdata, ProbClosures const &cls,
              ProbParm const &prob_parm) {
  const Real *prob_lo = geomdata.ProbLo();
  const Real *prob_hi = geomdata.ProbHi();
  const Real *dx = geomdata.CellSize();

  Real x = prob_lo[0] + (i + Real(0.5)) * dx[0];
  Real Pt, rhot, uxt;
  if (x < prob_hi[0] / 2) {
    Pt = prob_parm.p_l;
    rhot = prob_parm.rho_l;
    uxt = prob_parm.u_l;
  } else {
    Pt = prob_parm.p_r;
    rhot = prob_parm.rho_r;
    uxt = prob_parm.u_r;
  }
  state(i, j, k, cls.URHO) = rhot;
  state(i, j, k, cls.UMX) = rhot * uxt;
  state(i, j, k, cls.UMY) = Real(0.0);
  state(i, j, k, cls.UMZ) = Real(0.0);
  Real et = Pt / (cls.gamma - Real(1.0));
  state(i, j, k, cls.UET) = et + Real(0.5) * rhot * uxt * uxt;
}
```
