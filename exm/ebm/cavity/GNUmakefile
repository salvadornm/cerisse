TINY_PROFILE = FALSE

DEBUG    = FALSE
DIM      = 2
COMP     = gcc
CXXSTD   = c++20

USE_MPI = TRUE
USE_CUDA = FALSE
USE_OMP = FALSE
# Boundaries
USE_GPIBM=FALSE
USE_EB = TRUE

AMREX_NO_PROBINIT=TRUE

AMR_SOLVER = $(abspath ../../../)
AMREX_HOME = $(abspath ../../../lib/amrex)

CEXE_headers += prob.h
CEXE_sources += eb_geometry.cpp
include $(AMR_SOLVER)/src/Make.CNS

