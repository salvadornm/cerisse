TINY_PROFILE = FALSE

DEBUG   = FALSE
DIM     = 3
COMP     = gcc

USE_MPI = TRUE
USE_CUDA = FALSE
USE_OMP = FALSE
USE_GPIBM=FALSE

AMREX_NO_PROBINIT=TRUE
CUDA_ARCH=8.0

AMR_SOLVER ?= $(abspath ../../)
AMREX_HOME ?= $(abspath ../../lib/amrex)

include ./Make.package
include $(AMR_SOLVER)/src/Make.CNS
