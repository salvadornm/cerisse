# Running

This page explains details of the files needed to run the code, 
It is advisable to check the [Quickrun](quickrun.md) before reading this page.

## Building

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

See `amrex/Tools/GNUMake/` for a full list of supported compilers.

WARNING: Make sure that if gnu option is installed, g++ points to the right place. 
This happens often in Mac-OS, where g++ points to the native clang compiler (which is not supported).


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

`TINY_PROFILE` is the recommended way for profiling the code. When the executable is compiled with `TINY_PROFILE = TRUE`, it prints out something like this at the end of each execution, which gives you a good indication of how much computational time each portion of the code is taking up.

NOTE: You need to keep `PROFILE = FALSE` when setting `TINY_PROFILE = TRUE`.

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
They are passed to AMReX and PelePhysics

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
TEST = FALSE  # carry out additional tests, e.g. check if Fabs are properly initialised
DEBUG = FALSE # allow debugger to insert breakpoints, only support single-core
FSANITIZER = FALSE
THREAD_SANITIZER = FALSE
```

**CNS** options
  
```makefile
# CNS
USE_EB = FALSE
USE_FULL_PROB_POST_TIMESTEP = FALSE
USE_PROB_POST_COARSETIMESTEP = FALSE
USE_PROB_PARM_HOST = FALSE
USE_PMFDATA = FALSE
```
- `USE_EB`: compile with AMReX's embedded boundary method
- `USE_FULL_PROB_POST_TIMESTEP`: use `prob_post_timestep` or `CNS::full_prob_post_timestep`. The former is a kernel function that takes in variables like time, state, etc. and do something in a point-wise fashion, while the latter gives user more freedom, e.g., to access all members of the `CNS` class as well as perform complex multi-level operations.
- `USE_PROB_POST_COARSETIMESTEP`: use `CNS::prob_post_coarsetimestep` or not. It is similar to `CNS::full_prob_post_timestep`, but is called only after a coarse timestep.
- `USE_PROB_PARM_HOST`: define a `ProbParmHost` struct that behaves the same as the `ProbParm` but lives on host only. This is useful when user wants to store some large data and is not used on GPU, e.g. the initial conditions in the homogeneous isotropic turbulence (HIT) case.
- `USE_PMFDATA`: use `pmfdata` or not. This is a PelePhysics class that manages premixed flame (PMF) data generated by cantera. When this option is set to `TRUE`, users need to do `CNS::pmfdata.initialize()` in `amrex_probinit` to read the input .dat file, and then the `pmfdata` will be passed to `prob_initdata` and `bcnormal` for initialising the problem and setting boundary conditions (see PMF1D case for example).

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

You can find out all ready-made chemistry mechanisms at `PelePhysics/Support/Mechanism/Models`. If nothing suits you, consider building your own mechanism by converting the CHEMKIN or Cantera format files to PelePhysics files (see PelePhysics CETPR for more info).

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

Header file to define the problem. It defines 4 mandatory functions:

**`prob_initdata`** initialises the conservative state data:
```cpp
/**
 * \brief Initialise state data.
 * @param i          x position.
 * @param j          y position.
 * @param k          z position.
 * @param state      output state data.
 * @param geomdata   domain geometry data.
 * @param parm       Parm data defined in parm.H.
 * @param prob_parm  ProbParm data as defined in prob_parm.H and initialised in amrex_probinit.
 */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void prob_initdata(int i, int j, int k, amrex::Array4<amrex::Real> const& state,
                  amrex::GeometryData const& geomdata, Parm const& /*parm*/,
                  ProbParm const& prob_parm)
```

**`bcnormal`** sets the ghost cell states when the corresponding boundary condition is set to "Inflow / UserBC" in the input file (see [Input options](input-options.md)):
```cpp
/**
 * \brief Fill external boundary conditions for ghost cells.
 * @param x          ghost cell cooridinates.
 * @param s_int      flow state inside of the domain.
 * @param s_ext      flow state to be filled.
 * @param idir       direction (0: x, 1: y, 2: z).
 * @param sgn        high or low boundary (1: low, -1: high).
 * @param time       time.
 * @param geomdata   domain geometry data.
 * @param prob_parm  ProbParm data as defined in prob_parm.H and initialised in amrex_probinit.
 */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void bcnormal(const amrex::Real x[AMREX_SPACEDIM],
              const amrex::Real* /*s_int[LEN_STATE]*/,
              const amrex::Real* /*s_refl[LEN_STATE]*/, amrex::Real s_ext[LEN_STATE],
              const int idir, const int sgn, const amrex::Real /*time*/,
              amrex::GeometryData const& /*geomdata*/, ProbParm const& prob_parm)
```
NOTE: `prob_initdata` and `bcnormal` also take in a `pele::physics::PMF::PmfData::DataContainer const* pmf_data` if the compile-time option `USE_PMFDATA = TRUE`.

**`prob_post_restart`** runs after loading the data from checkpoint
```cpp
/**
 * \brief Modify state data and/or add turbulence to fields after restart.
 * @param i         x position.
 * @param j         y position.
 * @param k         z position.
 * @param state     state data.
 * @param geomdata  domain geometry data.
 * @param parm      Parm data defined in parm.H.
 * @param prob_parm ProbParm data as defined in prob_parm.H and initialised in amrex_probinit.
 */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE 
void prob_post_restart(int i, int j, int k, 
                       amrex::Array4<amrex::Real> const& state,
                       amrex::GeometryData const& geomdata, Parm const& parm, 
                       ProbParm const& prob_parm)
```

**`prob_post_timestep`** runs after each level timestep (not coarseTimestep). It can be used to record time averages:
```cpp
AMREX_GPU_DEVICE AMREX_FORCE_INLINE 
void prob_post_timestep(int i, int j, int k, 
                        amrex::Array4<amrex::Real> const& state,
                        amrex::GeometryData const& geomdata, Parm const& parm, 
                        ProbParm const& prob_parm)
```

Apart from standard refinement criteria provided in the input file, users can use **`prob_tag_error`** to define more complex criteria to flag regions for refinement 
```cpp
AMREX_GPU_DEVICE AMREX_FORCE_INLINE 
void prob_tag_error(int i, int j, int k, amrex::Array4<char> const& tagarr,
                    amrex::Array4<const amrex::Real> const& /*sarr*/, int level, char tagval,
                    const amrex::Real /*time*/, amrex::GeometryData const& geomdata,
                    Parm const& /*parm*/, ProbParm const& /*pp*/)
{
  if (/*some refinement criteria*/) {
    tagarr(i, j, k) = tagval; // refine
  }
}
```

NOTE: These functions must be defined. If they are not used in the case, just leave the function body empty. Omitting any will result in a compiler error.

If using auxiliary variables (e.g. time-averages), users must include the **`prob_get_aux_name`** to define the names of the auxiliary variables:
```cpp
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void prob_get_aux_name(amrex::Vector<std::string>& aux_name)
{
  aux_name.resize(NUM_AUX);
  aux_name[0] = "time_avg_u";
  aux_name[1] = "time_avg_v";
  aux_name[2] = "time_avg_w";
  // ...
}
```

## prob.cpp

We follow AMReX's convention of putting GPU/Device functions (those marked by `AMREX_GPU_DEVICE`) and declarations in `.H` files, and the CPU/Host functions in `.cpp` files.

The initialization workflow in the code is as follows: 

1. `amrex_probinit` initializes a ProbParm object `CNS::h_prob_parm`, which lives on CPU/Host
2. The data is copied over to GPU/Device ProbParm called `CNS::d_prob_parm`, at the end of `amrex_probinit`
3. `CNS::d_prob_parm` will be fed to `prob_initdata` to initialize the state data

The **`amrex_probinit`** function is defined in `prob.cpp`:

```cpp
extern "C" {
void amrex_probinit(const int* /*init*/, const int* /*name*/, const int* /*namelen*/,
                    const Real* /*problo*/, const Real* /*probhi*/)
{
  // Parse some options
  {
    amrex::ParmParse pp("prob");        // search for options with prob.XXX
    pp.get("T", CNS::h_prob_parm->T);   // pp.get produces error when option is not found 
    pp.query("p", CNS::h_prob_parm->p); // pp.quert does nothing when option is not found
  }

  // Do some init
  // ...

  amrex::Gpu::copy(Gpu::hostToDevice, CNS::h_prob_parm, CNS::h_prob_parm + 1,
                   CNS::d_prob_parm);
}
}
```

It usually contains the `amrex::ParmParse` to parse input file options.

NOTE: Some examples use `Gpu::copyAsync()` followed by `Gpu::streamSynchronize()`. This is the same as `Gpu::copy()`, which implicitly runs `Gpu::streamSynchronize()`. If one wants to copy more than one object to/from GPU, the `Gpu::copyAsync()` followed by `Gpu::streamSynchronize()` method is preferred as only synchronization will be called.

The `fill_ext_src` function is called if the `cns.do_ext_src` option is true. It can be used to add gravity or ignition.

```cpp
void CNS::fill_ext_src(int i, int j, int k, amrex::Real time,
                       amrex::GeometryData const& geomdata,
                       amrex::Array4<const amrex::Real> const& state,
                       amrex::Array4<amrex::Real> const& ext_src,
                       Parm const& /*parm*/, ProbParm const& pp)
{
  // Add some source term here, e.g.
  Real uz = state(i, j, k, UMZ) / state(i, j, k, URHO);
  ext_src(i, j, k, UMZ) += g;
  ext_src(i, j, k, UEDEN) += uz * g;
}
```

DANGER: The `ext_src` array is already filled with hydrodynamic and viscous sources. Make sure you use `+=` instead of `=` when adding external sources. This design is intended to give developers more control, e.g. can add source terms as a function of the time derivative of the states. (Maybe we should change this in the future?)

## prob_param.H

This header file defines the `ProbParm` class which holds all the problem-specific data for initialization, boundary conditions, etc. As mentioned before, it is accessible to both Host and Device, so all data must be GPU-safe. For example, use `amrex::GpuArray` instead of pointer-type `amrex::Real []` for storing arries.

```cpp
struct ProbParm
{
  amrex::Real p = 1.0;
  amrex::Real u = 0.0;
  amrex::Real rho = 1.0;
  amrex::Real rhoe;
  amrex::Real T;
  amrex::GpuArray<amrex::Real, NUM_SPECIES> massfrac = {0.0};
};
```

Note that the `= {0.0}` syntax initialises every entry of the array to zero. It only works for zeros. See the [C++ reference](https://en.cppreference.com/w/c/language/array_initialization) for more information.


## input file

The input file contains runtime (post-compile time) parameters. See [Input options](input-options.md) for full information. 


## Embedded boundary geometries

Cerisse supports AMReX's cut-cell embedded boundary method for geometries. Users can use the built-in AMReX shapes (e.g. "box", "cylinder", "plane", "sphere", "torus"), or Cerisse's "triangles", or read the geometry from an stl file. 

Alternatively, users can also build their geometries using the combination of the built-in AMReX shapes, e.g. in the Michigan Scramjet case:

First, declare the geometry in `prob.H`. Here, `Scamjet` is the name of the geometry, replace that with your own name.
```cpp
#include "custom_geometry.H"
class Scramjet : public CustomGeometry::Register<Scramjet>
{
public:
  static const std::string identifier() { return "Scramjet"; }

  void build(const Geometry& geom, const int max_coarsening_level) override;
};
```
TIP: The `CustomGeometry` class is an application of the PelePhysics `Factory` class. It is a neat way to create options for different solution methods. See `PelePhysics/Source/Factory.H` for more information.

Then, define the `build` method in the `prob.cpp` file. Note that the function definition cannot be in `prob.H` or else it will be defined multiple times, which will result in a compiler error.
```cpp
void Scramjet::build(const Geometry& geom, const int max_coarsening_level)
{
  // Step by step, we first create the part before the cavity, which 
  // is a box with a cylinder hole as the injector. We use DifferenceIF 
  // to subtract the cylinder from the box
  auto box = EB2::BoxIF({AMREX_D_DECL(-50., -10., -10.)}, {AMREX_D_DECL(4.45, 0.0, 10.)}, false);
  auto injector = EB2::CylinderIF(0.1245, 2.0, 1, {AMREX_D_DECL(0.0, -0.635, 0.0)}, false);
  auto box_with_inj = EB2::DifferenceIF<EB2::BoxIF, EB2::CylinderIF>(box, injector);

  // Then, the part after the cavity is a 4 degree wedge. We create 
  // 3 planes and intersect them together to form a triangle (the part 
  // extended outside the domain will be obmitted so it will give us 
  // the wedge). Notice that the wall normal must point into the solid.
  auto rear_wall = EB2::PlaneIF({AMREX_D_DECL(9.525, 0.0, 0.0)}, {AMREX_D_DECL(1.0, 0.0, 0.0)});
  auto floor_wall = EB2::PlaneIF({AMREX_D_DECL(0.0, -10.0, 0.0)}, {AMREX_D_DECL(0.0, 1.0, 0.0)});
  auto inclined_wall = EB2::PlaneIF({AMREX_D_DECL(9.525, 0.0, 0.0)}, {AMREX_D_DECL(-sin(4.0 / 180.0 * M_PI), -cos(4.0 / 180.0 * M_PI), 0.0)});
  auto triangle = EB2::IntersectionIF<EB2::PlaneIF, EB2::PlaneIF, EB2::PlaneIF>(rear_wall, floor_wall, inclined_wall);
    
  // Finally, combine the two parts together using UnionIF
  auto all_objs = EB2::UnionIF<EB2::DifferenceIF<EB2::BoxIF, EB2::CylinderIF>, 
                  EB2::IntersectionIF<EB2::PlaneIF, EB2::PlaneIF, EB2::PlaneIF>>(box_with_inj, triangle);
  auto gshop = EB2::makeShop(all_objs);
  EB2::Build(gshop, geom, max_coarsening_level, max_coarsening_level, 6, true);
}
```

The basic geometries (implicit functions or IF) supported by AMReX are:
- **Plane** - needs a point (plane_point) and normal (plane_normal).
- **Sphere** - needs center (sphere_center), radius (sphere_radius) and fluid inside/outside flag (sphere_has_fluid_inside).
- **Cylinder** - needs center (cylinder_center), radius (cylinder_radius), height (cylinder_height), direction (cylinder_direction) and fluid inside/outside flag (cylinder_has_fluid_inside).
- **Box** - needs the lower corner (box_lo), upper corner (box_hi) and fluid inside/outside flag (box_has_fluid_inside).

And the transformation handles in AMReX are:
- **Intersection** - find the common region between implicit functions
- **Union** - find the union of implicit functions
- **Complement** - invert an implicit function, i.e. make fluid that is inside to outside
- **Translation** - translate an implicit function
- **Lathe** - creates a 3D implicit function from a 2D function by revolving about the z-axis
- **Extrusion** - creates a 3D implicit function from a 2D function by translating along the z-axis