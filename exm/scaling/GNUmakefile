TINY_PROFILE 	= FALSE
DEBUG   	= FALSE
PRECISION       = DOUBLE 
DIM     	= 3
COMP 		= gcc
CXXSTD    	= c++20
USE_CUDA 	= FALSE 
USE_MPI 	= TRUE
USE_OMP 	= FALSE
USE_GPIBM	= FALSE
AMREX_NO_PROBINIT=TRUE
CUDA_ARCH=8.0

AMR_SOLVER = $(abspath ../../)
AMREX_HOME = $(abspath ../../lib/amrex)

CEXE_headers += prob.h
include $(AMR_SOLVER)/src/Make.CNS


# Temporarily kill heavy optimization even when DEBUG=FALSE
#XTRA_CXXFLAGS += -O1
#XTRA_FFLAGS   += -O1

# 2) Sanitizers to catch memory / UB issues
#XTRA_CXXFLAGS += -fsanitize=address,undefined -fno-omit-frame-pointer
#XTRA_LDFLAGS  += -fsanitize=address,undefined
#XTRA_CXXFLAGS += -fsanitize=undefined -fno-omit-frame-pointer
#XTRA_LDFLAGS  += -fsanitize=undefined

# Try disabling strict aliasing in optimized build  for Production runs
XTRA_CXXFLAGS += -fno-strict-aliasing
XTRA_FFLAGS   += -fno-strict-aliasing



