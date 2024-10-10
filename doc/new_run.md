# CODE DISCUSSION DOCUMENT


The code is based on two premises:

### Physical Problem Closures

creating the "closures" of the problems, that means which model of transport 
properties is to be used. This depend on the problem and is seletced in `prob.h`

For example:

```cpp
typedef closures_dt<indicies_t, visc_suth_t, cond_suth_t,
                    calorifically_perfect_gas_t<indicies_t>> ProbClosures;
```

This selects a problem that uses the structure `indicies_t`, which tell position
of the variables and how many variables to solve. The above also selects
Sutherland viscosity model  ```visc_suth_t```
All this is wraped in the class ```ProbClosures``


### Equations to Solve and how

Define the RHS of the problem, that includes which equation to solve and which nuemrivcal scheme to use. In general 

$$
\frac{\partial U}{\partial t} = \mbox{RHS}(U)
$$

where the RHS includes, inviscid (Euler) terms viscous (for Navier-Stokes) and source terms. For example:

```cpp
typedef rhs_dt<skew_t<true,false, 4, ProbClosures>, no_diffusive_t, no_source_t > ProbRHS;
```
This will solve
$$
\frac{\partial U}{\partial t} = F(U)
$$
corresponding to the Euelr equations, with a 4th order skew-symmetric numerical scheme.
Template `skew_t`,`diffusive_t` may have  arguments (like the order of the scheme).

**NEW** Parameters can be passed by a small structure `methodparm_t`,
so all paraeters  can be changed directly in `prob.h`, 
such as the order of the scheme, disspation parameters, 
extra damping, sensor variables, etc. The above example would be

```cpp
struct methodparm_t {

  public:

  static constexpr int  order = 2;              // order numerical scheme viscous
  static constexpr Real conductivity = 0.0262;  // conductivity (for constant value)
  static constexpr Real viscosity   = 1.85e-5;  // viscosity    (for constant value)
  
};
typedef closures_dt<indicies_t, visc_const_t<methodparm_t>, cond_const_t<methodparm_t>,
                    calorifically_perfect_gas_t<indicies_t>> ProbClosures;

typedef rhs_dt<no_euler_t, diffusiveheat_t<methodparm_t, ProbClosures>, no_source_t > ProbRHS;

```

This is used to set manually the conductivity and viscosity to a desired value ensure 
and accuarcy of derivatives in `diffusiveheat` 


Variables


Conservative

rhok,


Primitive

Storing but losing




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

This is done in the file **prob.h**.
In the top of the file, the following line are always the same:


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

using the use AMREX geometry and call to Closures and RHS 
(link here to pages).

This is follow by namespace ```PROB {```
and the definition of problem parameters:

```cpp
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

that may be useful.
The names follow in the **cons_var_names** array

```cpp
inline Vector<std::string> cons_vars_names={"Xmom","Ymom","Zmom","Energy","Energy"};
```

The type of variables, keep as it is, scalar ser 0 and vectors are given by their components


```cpp
inline Vector<int> cons_vars_type={1,2,3,0,0};
```

The closures and right-hand-side of equations

```cpp
  typedef closures_dt<indicies_t, visc_suth_t, cond_suth_t,
                    calorifically_perfect_gas_t<indicies_t>>
    ProbClosures;
  typedef rhs_dt<riemann_t<false, ProbClosures>, no_diffusive_t,no_source_t>
    ProbRHS;
```

### Closures


| Closure                      | Type          | Use | Description                                                  |
| --------------------------- | ------------- |:-------:| ------------------------------------------------------------ |
| ```visc_suth```             |            |   often      | Sutherland Viscosity                           |
| ```cond_suth```             |            |   often      | Sutherland Conductivty                         |
| ```calorifically_perfect_gas_t```             |            |   often      | perfect gas                      |


### RHS

This call indicates the different type of numerical scheme for the advection part

```cpp
  typedef rhs_dt<keep_euler_t<false, false, 4, ProbClosures>,no_diffusive_t,
               no_source_t> ProbRHS;
```

| RHS                     | Options          | Use | Description                                                  |
| --------------------------- | ------------- |:-------:| ------------------------------------------------------------ |
| ```riemann_t```             |            |   often      | Riemann Solver   MUSCL              |
| ```keep_euler_t```             |  AD,IB, order           |   often      | KEEP Scheme                 |
| ```centraldif_t```             |  AD,IB, order           |   often      | Central Scheme                 |
| ```skew_t```             |  AD, IB, order           |   often      | Skew-symmetric Scheme                 |
| ```no_diffusive_t```        |            |   often      | No diffusive part (Euler)                |
| ```no_source_t```           |            |   often      | No source term    |



### Inputs


Data missed from input file

```cpp
void inline inputs() {
  ParmParse pp;

  pp.add("cns.order_rk", 3);   // -2, 1, 2 or 3"
  pp.add("cns.stages_rk", 3);  // 1, 2 or 3
}
```

## Problem definition

The following lines define the initial conditions

```cpp
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
