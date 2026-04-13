# IBM test: CGAL CPU 2D (single circle)
DEBUG     = FALSE
DIM       = 2
COMP      = gcc
CXXSTD    = c++20
PRECISION = DOUBLE
TINY_PROFILE = TRUE
CLIP_MINTEMP = FALSE

USE_MPI  = TRUE
USE_CUDA = FALSE
USE_OMP  = FALSE
USE_GPIBM = TRUE
USE_CGAL  = TRUE

GMP_HOME  = /usr
MPFR_HOME = /usr

AMREX_NO_PROBINIT = TRUE
CUDA_ARCH = 8.0

AMR_SOLVER = $(abspath ../../../../)
AMREX_HOME = $(AMR_SOLVER)/lib/amrex

CEXE_headers += prob.h
VPATH_LOCATIONS += $(abspath ../)
INCLUDE_LOCATIONS += $(abspath ../)
include $(AMR_SOLVER)/src/Make.CNS
