
DEBUG   = FALSE
DIM     = 3
COMP     = gcc
CXXSTD    = c++20
PRECISION = DOUBLE

TINY_PROFILE = FALSE

CLIP_MINTEMP = TRUE

USE_MPI = TRUE
USE_CUDA = FALSE
USE_OMP = FALSE

# Boundaries
USE_GPIBM=FALSE
USE_EB = TRUE

AMREX_NO_PROBINIT=TRUE
CUDA_ARCH=8.0

AMR_SOLVER = $(abspath ../../../)
AMREX_HOME = $(abspath ../../../lib/amrex)

CEXE_headers += prob.h
CEXE_sources += eb_geometry.cpp
include $(AMR_SOLVER)/src/Make.CNS
