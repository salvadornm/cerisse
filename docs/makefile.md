---
icon: laptop-code
cover: >-
  https://images.unsplash.com/photo-1577401239170-897942555fb3?crop=entropy&cs=srgb&fm=jpg&ixid=M3wxOTcwMjR8MHwxfHNlYXJjaHw2fHxwcm9ibGVtfGVufDB8fHx8MTczMDk3NjczN3ww&ixlib=rb-4.0.3&q=85
coverY: 0
---

# Makefile

This file sets general options, it is usually modified once and that are required to compile the code providing options to AMReX, PelePhysics and cerisse itself. As `prob.h` is best to re-use files from examples The file options can divided into groups

## AMRex/Cerisse options

These lines are essential for configuring AMReX and setting the simulation dimensions, which is important for optimizing space and enabling efficient 1D or 2D simulations. They also specify the compiler, precision, and the use of parallelization frameworks like MPI, OpenMP, and GPU acceleration with CUDA. Enabling debugging will compile with additional checks (significantly slowing down execution), while profiling options provide a summary of exclusive and inclusive function times, as well as the minimum and maximum time spent in each routine across processes. For more infromation, see [AMReX Profiling](https://amrex-codes.github.io/amrex/docs\_html/AMReX\_Profiling\_Tools.html)

```makefile
# AMReX
DIM = 1             # dimension of the problem (1/2/3)
COMP = gnu          # compiler (gnu/intel/..)
PRECISION = DOUBLE  # floating-point precision (it has not been tested in SINGLE)

# Performance
USE_MPI = TRUE
USE_OMP = FALSE
USE_CUDA = FALSE
CUDA_ARCH= 8.0

# Debugging
DEBUG = FALSE

# Profiling
TINY_PROFILE = FALSE
```

The below options related to model that involve compiling different files.

```makefile
USE_GPIBM = FALSE           # use of Immersed Boundaries for solid boundaries
USE_EB = FALSE              # use of Emmbedded Boundaries for solid boundaries
USE_PELEPHYSICS = TRUE      # PelePhysics is used for chemistry/thermodynamics/tarnsport
```



{% hint style="info" %}
Embedded Boundaries are **AMReX-native**, whiel Immersed Boundaries are **Cerisse**. Both are incompatible (check Boundary Conditions)
{% endhint %}

## PelePhysics

If PelePhysics selected, the GNU Makefile specifies the thermodynamics, transport, and chemistry mechanisms. It uses the [PelePhysics](https://pelephysics.readthedocs.io/en/latest/index.html) style.

{% hint style="warning" %}
The first time a PelePhysics simulation is running do `make TPL` to download and compile third-party libraries
{% endhint %}

```makefile
# PelePhysics
EOS_MODEL := FUEGO
TRANSPORT_MODEL := SIMPLE
CHEMISTRY_MODEL := BurkeDryer
```

There are plenty (around 30) of chemical mechanims stored in `./lib/PelePhysics/Support/Mechanism/Models` and more can be converted from [Cantera](https://cantera.org) format.

## Paths and global files

These variables specify the paths for AMReX and third-party libraries. The paths should either be local or, for more flexibility, use cleverly defined global paths.

```makefile
AMR_SOLVER ?= $(abspath ../../)
AMREX_HOME ?= $(abspath ../../lib/amrex)
PELE_PHYSICS_HOME = $(abspath ../../lib/PelePhysics)
```

The lines below add the user-specified `prob.h` (or different header file), as well as the path for the overall `Makefile`

```makefile
CEXE_headers += prob.h
include $(AMR_SOLVER)/src/Make.CNS
```
