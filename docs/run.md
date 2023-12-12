# Running

This page explains details of the files needed to run the code, 
It is advisable to check the [Quickrun](quickrun.md)  before reading this page.


### GNUmakefile 

This files sets general options, usually modified once and that are required to compile the code.
It is divided into sections

**AMREX** options
such as dimension of the problem (1/2/3), compiler (gnu/intel/nvc/..) and precision (single/double).

WARNING: Make sure that if gnu option is installed, g++ points to the right place. 
This happens often in Mac-OS, where g++ points to native clang compiler (which is not supported).


```
# AMReX
DIM = 1
COMP = gnu
PRECISION = DOUBLE
```

**Profiling** options


**Performance** options related to parallelization, using MPI/OMP/CUDA for GPU, etc.
This are passed to AMReX and PelePhysics

```
USE_MPI = TRUE
USE_OMP = FALSE
USE_CUDA = FALSE
USE_HIP = FALSE
USE_SYCL = FALSE
```

**Debugging** options

```
# Debugging
DEBUG = FALSE
FSANITIZER = FALSE
THREAD_SANITIZER = FALSE
```

**PelePhysisc** options, related to thermodynamics models, transport and chemistry

```
# PelePhysics
Eos_Model := GammaLaw
Transport_Model := Constant
Chemistry_Model := Null
```

variable |  meaning  | options
:----------- |:-------------:| -----------:
Eos_Model       | Equation of State        | GammaLaw/Fuego
Transport_Model     | Transport Model        | Constant/Simple
Chemistry_Model        | Chemistry Mechanism    | Null/grimech30/LiDryer



**Makefile** options, related to files to add to makefile

```
# GNU Make
include ./Make.package
include ../Make.CNS%
```

### Make.package

This file just add the names of the "problem" files to compile to headers and source files

```
CEXE_headers += prob.H prob_parm.H
CEXE_sources += prob.cpp
```


### input file

Frequent changes, in some files is divided into sections.
Comments can be added by `#`


```
max_step = 100
stop_time = 0.2
# ================ COMPUTATIONAL DOMAIN ================ #
geometry.is_periodic = 0 0 0
geometry.coord_sys   = 0  # 0 => cart, 1 => RZ  2=>spherical
geometry.prob_lo     = 0.    0.     0.
geometry.prob_hi     = 1.    0.25   0.25
amr.n_cell           = 32    8      8
```

Input file options

variable | type | meaning
:----------- |:-------------:| -----------:
max_step       | integer        | number of steps of simulations
stop_time       | real        | physical time to stop simulation
geometry.is_periodic        | 3*integer       | indicates periodic direction 
geometry.coord_sys  | integer |   cartesian cylindrical spherical
geometry.prob_lo   |  3*real    | domain lowersize x, y, z
geometry.prob_hi   |  3*real    | domain uppersize x, y, z
amr.n_cell         |  3*integer |   mesh in x,y,z (coarsest level) 

Self explanatory

```
# ================ BOUNDARY CONDITIONS ================= #
# 0 = Interior             3 = Symmetry
# 1 = Inflow / UserBC      4 = SlipWall =3
# 2 = Outflow (FOExtrap)   5 = NoSlipWall (adiabatic)
cns.lo_bc = 0 0 0
cns.hi_bc = 0 0 0
```

variable | type | meaning
:----------- |:-------------:| -----------:
cns.lo_bc       | 3*integer        | bc flags lower boundaries in x,y,z
cns.hi_bc       | 3*integer        | bc flags upper boundaries in x,y,z


```
# ==================== WHAT NUMERICS =================== #
cns.cfl = 0.1
cns.dt_cutoff = 5.e-20
cns.recon_scheme = 6  # 1: Godunov,  2: MUSCL,   3: WENO-Z3, 
                      # 4: WENO-JS5, 5: WENO-Z5, 6: TENO-5
```
Godunov is first order recontruction, MUSCL is second order, while the
number represent the formal order of the reconstruction scheme.


variable | type | meaning
:----------- |:-------------:| -----------:
cns.cfl       | real        | Acoustic Courant Number 
cns.dt_cutoff      | real        | bc flags upper boundaries in x,y,z
cns.recon_scheme   |  integer   | numerical scheme 
 


### prob.h

Header file to define problem

```
#ifndef CNS_PROB_H
#define CNS_PROB_H

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_REAL.H>

#include "CNS.H"
#include "prob_parm.H"

/**
 * \brief Initialise state data.
 *
 * @param i         x position.
 * @param j         y position.
 * @param k         z position.
 * @param state     state data.
 * @param geomdata  domain geometry data.
 * @param parm      Parm data defined in parm.H.
 * @param prob_parm ProbParm data as defined in prob_parm.H and initialised in
 * amrex_probinit.
 */
```

### prob.cpp


Start of file

```
#include "prob.H"

using namespace amrex;

extern "C" {
void amrex_probinit(const int* /*init*/, const int* /*name*/, const int* /*namelen*/,
                    const Real* /*problo*/, const Real* /*probhi*/)
{
}
}
```

### prob_param.H

An additional parameter file can be used, for example


```
#ifndef CNS_PROB_PARM_H
#define CNS_PROB_PARM_H

#include <AMReX_GpuMemory.H>
#include <AMReX_REAL.H>

using namespace amrex::literals;

struct ProbParm
{
  amrex::Real p_l = 1.0;   // left pressure
  amrex::Real u_l = 0.0;   // left velocity
  amrex::Real rho_l = 1.0; // left density
  amrex::Real rhoe_l;
  amrex::Real T_l;

  amrex::Real p_r = 0.1;     // right pressure
  amrex::Real u_r = 0.0;     // right velocity
  amrex::Real rho_r = 0.125; // right density
  amrex::Real rhoe_r;
  amrex::Real T_r;
};

#endif
```

## Geometry EB

In **input** file

```
##---- GEOMETRY -----
#eb2.geom_type = all_regular

eb2.geom_type = cylinder
eb2.cylinder_direction = 2
eb2.cylinder_radius = 0.25
eb2.cylinder_center = 1.0 2.0 0.0
eb2.cylinder_has_fluid_inside = 0

#eb2.sphere_radius = 0.5
#eb2.sphere_center = 2.0 2.0  2.0
#eb2.sphere_has_fluid_inside = 0

#eb2.geom_type = box
#eb2.box_lo = 0.75  1.75  0.0
#eb2.box_hi = 1.25  2.25  0.0
#eb2.box_has_fluid_inside = 0
```

There are several basic geometries that are available in AMReX that can be easily specified in the input file:
- Box
- Plane
- Cylinder
- Sphere

Check [PeleC](https://amrex-combustion.github.io/PeleC/geometry/EB.html) 
page about EB geometry 



## Options




