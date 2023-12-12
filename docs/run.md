# Running

This page explains details of the files needed to run the code, 
It is advisable to check the [Quickrun](quickrun.md)  before reading this page.

## Building

There are two ways to compile the code, using **GNUMake** or **CMake**. The former is recommended because it is simple and more comprehensible, while the latter is mainly used for automated testing with CTest.

### Required files

Each case must contain the following 5 files to get compiled and run

* **`GNUmakefile`** - sets the compile-time options
* **`prob_parm.H`** - define the `ProbParm` struct, which contains data for initialization or boundary conditions
* **`prob.H`** - define functions for initialization (`prob_initdata`), boundary conditions (`bcnormal`), etc.
* **`prob.cpp`** - define function to initialize AMReX (`amrex_probinit`)
* **`inputs`** - contain runtime options for the cerisse executable. Some may contain the `.ini` extension, but they are the same thing

There may also be a few optional files

* **`Make.package`** - link the `.H` and `.cpp` files to the compiler. It can be absorbed into `GNUmakefile`
* **`CMakeLists.txt`** - for building with CMake

We will now go through these files one by one.

## GNUmakefile

This files sets general options, usually modified once and that are required to compile the code.
It is divided into sections

**AMReX** options

```makefile
# AMReX
DIM = 1             # dimension of the problem (1/2/3)
COMP = gnu          # compiler (gnu/intel/..)
PRECISION = DOUBLE  # floating-point precision (it has not been tested in SINGLE)
```

See `amrex/Tools/GNUMake/` for a full list of supported compilers.

WARNING: Make sure that if gnu option is installed, g++ points to the right place. 
This happens often in Mac-OS, where g++ points to native clang compiler (which is not supported).


**Profiling** options

```makefile
# Profiling
PROFILE = FALSE 
TINY_PROFILE = FALSE  # use amrex light-weight profiler (recommended)
COMM_PROFILE = FALSE
TRACE_PROFILE = FALSE
MEM_PROFILE = FALSE
USE_GPROF = FALSE
```

`TINY_PROFILE` is the recommended way for profiling the code. When the executable is compiled with `TINY_PROFILE = TRUE`, it prints out something like this at the end of each execution, which gives you a good indication of what how much of computational time each portion of the code is taking up.

NOTE: You need to keep `PROFILE = FALSE` when `TINY_PROFILE = TRUE`.

```bash
TinyProfiler total time across processes [min...avg...max]: 2.383 ... 2.383 ... 2.383

--------------------------------------------------------------------------------------------
Name                                         NCalls  Excl. Min  Excl. Avg  Excl. Max   Max %
--------------------------------------------------------------------------------------------
Pele::ReactorRK64::react()                      172      1.024      1.195      1.322  55.47%
FabArray::ParallelCopy_finish()                2050     0.1004      0.282     0.5157  21.64%
FillBoundary_finish()                           224     0.2513     0.3731     0.4651  19.52%
...
--------------------------------------------------------------------------------------------

--------------------------------------------------------------------------------------------
Name                                         NCalls  Incl. Min  Incl. Avg  Incl. Max   Max %
--------------------------------------------------------------------------------------------
main()                                            1      2.383      2.383      2.383 100.00%
Amr::coarseTimeStep()                            10      2.288      2.288      2.288  96.04%
Amr::timeStep()                                 150      2.285      2.286      2.286  95.94%
...
--------------------------------------------------------------------------------------------

Pinned Memory Usage:
---------------------------------------------------------------------------------------------------------------------
Name                            Nalloc  Nfree  AvgMem min  AvgMem avg  AvgMem max  MaxMem min  MaxMem avg  MaxMem max
---------------------------------------------------------------------------------------------------------------------
The_Pinned_Arena::Initialize()       8      8    9346   B      17 KiB      28 KiB    8192 KiB    8192 KiB    8192 KiB
---------------------------------------------------------------------------------------------------------------------
```

**Performance** options related to parallelization, using MPI/OMP/CUDA for GPU, etc.
This are passed to AMReX and PelePhysics

```makefile
# Performance
USE_MPI = TRUE
USE_OMP = FALSE
USE_CUDA = FALSE
USE_HIP = FALSE
USE_SYCL = FALSE
```

**Debugging** options

```makefile
# Debugging
DEBUG = FALSE
FSANITIZER = FALSE
THREAD_SANITIZER = FALSE
```

**PelePhysisc** options, related to thermodynamics models, transport and chemistry

```makefile
# PelePhysics
Eos_Model := GammaLaw
Transport_Model := Constant
Chemistry_Model := Null
```

Variable        | Meaning             | Options
:-------------- |:--------------------| :-------------------------------
Eos_Model       | Equation of State   | GammaLaw / Fuego / SRK
Transport_Model | Transport Model     | Constant / Simple
Chemistry_Model | Chemistry Mechanism | Null / grimech30 / LiDryer / ...

You can find out all ready-made chemistry mechanisms at `PelePhysics/Support/Mechanism/Models`. If nothing suits, consider building your own mechanism by converting the CHEMKIN or Cantera format files to PelePhysics files.

**Makefile** options, related to files to add to makefile

```makefile
# GNU Make
include ./Make.package  # insert the Make.package contents here
include ../Make.CNS     # include rules to build the executable
```

You may choose to remove the `Make.package` file and manually insert all its contents into `GNUmakefile`, like this:

```makefile
# GNU Make
CEXE_headers += prob.H prob_parm.H  # These are in the Make.package
CEXE_sources += prob.cpp            #
include ../Make.CNS     # include rules to build the executable
```

## prob.H

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

## prob.cpp


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

## prob_param.H

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

## input file

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




