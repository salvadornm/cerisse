# FSI: Moving diamond wedge airfoil — CPU build
DEBUG     = FALSE
DIM       = 2
COMP      = gcc
CXXSTD    = c++20
PRECISION = DOUBLE
TINY_PROFILE = TRUE
CLIP_MINTEMP = TRUE

USE_MPI  = TRUE
USE_CUDA = FALSE
USE_OMP  = FALSE

USE_GPIBM = TRUE
USE_FSI   = TRUE

AMREX_NO_PROBINIT = TRUE
CUDA_ARCH = 8.0

AMR_SOLVER = $(abspath ../../../)
AMREX_HOME = $(AMR_SOLVER)/lib/amrex

CEXE_headers += prob.h
include $(AMR_SOLVER)/src/Make.CNS
