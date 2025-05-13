DEBUG   = FALSE
DIM     = 3
COMP     = gcc
CXXSTD    = c++20
PRECISION = DOUBLE

TINY_PROFILE = FALSE

# cerisse options
CLIP_MINTEMP = FALSE

# parallelization
USE_MPI = TRUE
USE_CUDA = FALSE
USE_OMP = FALSE

# PelePhysics
USE_PELEPHYSICS = TRUE
EOS_MODEL := FUEGO
TRANSPORT_MODEL := SIMPLE
CHEMISTRY_MODEL := MartianAtmos
USE_LOCALCHEM = TRUE
LOCALCHEM_PATH = $(abspath ./)

# Boundaries
USE_GPIBM=TRUE

AMREX_NO_PROBINIT=TRUE
CUDA_ARCH=8.0

# Paths
AMR_SOLVER = $(abspath ../../../)
AMREX_HOME = $(abspath ../../../lib/amrex)
PELE_PHYSICS_HOME = $(abspath ../../../lib/PelePhysics)


CEXE_headers += prob.h
include $(AMR_SOLVER)/src/Make.CNS
