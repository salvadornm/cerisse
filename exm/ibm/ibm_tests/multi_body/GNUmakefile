# IBM test: Multi-body (7 cylinders, BVH CPU 2D)
DEBUG     = FALSE
DIM       = 2
COMP      = gcc
CXXSTD    = c++20
PRECISION = DOUBLE
TINY_PROFILE = FALSE
# Enable T_min floor: essential for Mach-4 shock/body interaction to avoid
# negative internal-energy overshoot in WENO/Roe reconstruction near bodies.
CLIP_MINTEMP = TRUE

USE_MPI  = TRUE
USE_CUDA = FALSE
USE_OMP  = FALSE
USE_GPIBM = TRUE

AMREX_NO_PROBINIT = TRUE
CUDA_ARCH = 8.0

AMR_SOLVER = $(abspath ../../../../)
AMREX_HOME = $(AMR_SOLVER)/lib/amrex

CEXE_headers += prob.h
include $(AMR_SOLVER)/src/Make.CNS
