TINY_PROFILE = FALSE
DEBUG = FALSE
DIM=3
USE_MPI = TRUE
USE_CUDA= FALSE
USE_OMP = FALSE
USE_GPIBM=TRUE
AMREX_NO_PROBINIT=TRUE
CUDA_ARCH=8.0
COMP     = gcc

AMR_SOLVER = $(abspath ../../)
AMREX_HOME = $(abspath ../../lib/amrex)

CEXE_headers += prob.h
include $(AMR_SOLVER)/src/Make.CNS
