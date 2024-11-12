TINY_PROFILE 	= FALSE
DEBUG   	= FALSE
PRECISION       = DOUBLE 
DIM     	= 3
COMP 		  = gcc
CXXSTD    = c++20
USE_CUDA 	= FALSE 
USE_MPI 	= TRUE
USE_OMP 	= FALSE
USE_GPIBM	= FALSE
AMREX_NO_PROBINIT=TRUE
CUDA_ARCH=8.0

AMR_SOLVER = $(abspath ../../../)
AMREX_HOME = $(abspath ../../../lib/amrex)

CEXE_headers += prob.h
include $(AMR_SOLVER)/src/Make.CNS
