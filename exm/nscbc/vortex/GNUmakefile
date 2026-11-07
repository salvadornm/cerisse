TINY_PROFILE = FALSE

DEBUG    = FALSE
DIM      = 2
COMP     = gcc
CXXSTD   = c++20

USE_MPI = TRUE
USE_CUDA = FALSE
USE_OMP = FALSE
USE_GPIBM=FALSE

AMREX_NO_PROBINIT=TRUE

AMR_SOLVER = $(abspath ../../../)
AMREX_HOME = $(abspath ../../../lib/amrex)

include ./Make.package
include $(AMR_SOLVER)/src/Make.CNS
