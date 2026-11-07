# IBM test: BVH GPU 3D (single sphere)
DEBUG     = FALSE
DIM       = 3
COMP      = gcc
CXXSTD    = c++20
PRECISION = DOUBLE
TINY_PROFILE = TRUE
CLIP_MINTEMP = TRUE

USE_MPI  = TRUE
USE_CUDA = TRUE
USE_OMP  = FALSE
USE_GPIBM = TRUE

AMREX_NO_PROBINIT = TRUE
CUDA_ARCH = 8.0

AMR_SOLVER = $(abspath ../../../../)
AMREX_HOME = $(AMR_SOLVER)/lib/amrex

CEXE_headers += prob.h
VPATH_LOCATIONS += $(abspath ../)
INCLUDE_LOCATIONS += $(abspath ../)
include $(AMR_SOLVER)/src/Make.CNS
