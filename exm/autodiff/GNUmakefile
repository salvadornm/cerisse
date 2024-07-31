TINY_PROFILE 	= FALSE
DEBUG   	= FALSE
PRECISION = DOUBLE # DOUBLE or FLOAT
DIM     	= 3
COMP 		  = llvm #GNU # we need to use llvm for clad
CXXSTD    = c++20
USE_CUDA 	= TRUE	
USE_MPI 	= FALSE
USE_OMP 	= FALSE
USE_GPIBM	= FALSE
AMREX_NO_PROBINIT=TRUE
CUDA_ARCH=8.0

AMR_SOLVER = $(abspath ../../)
AMREX_HOME = $(abspath ../../lib/amrex)

INCLUDE_LOCATIONS += $(AMR_SOLVER)/lib/install/clad/include 
NVCC_CCBIN=/usr/bin/clang++

CEXE_headers += prob.h
include $(AMR_SOLVER)/src/Make.CNS