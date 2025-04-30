# AMReX
DIM = 2
COMP = gnu
CXXSTD    = c++20
PRECISION = DOUBLE

# Profiling
TINY_PROFILE = FALSE

# Performance
USE_MPI = TRUE
USE_OMP = FALSE
# USE_CUDA = TRUE
# CUDA_ARCH= 8.0

# Debugging
DEBUG = FALSE

# PelePhysics
EOS_MODEL := FUEGO
TRANSPORT_MODEL := SIMPLE
CHEMISTRY_MODEL := LiDryer

# Utilities
USE_UTILITIES = TRUE

# CNS
AMREX_NO_PROBINIT = TRUE
USE_GPIBM = FALSE
USE_EB = FALSE
USE_PELEPHYSICS = TRUE
#USE_SUNDIALS_PP = TRUE #for CVODE on GPU -- NOT needed https://github.com/AMReX-Combustion/PelePhysics/pull/194
#PELE_USE_KLU = FALSE #for direct solve

# GNU Make
AMR_SOLVER = $(abspath ../../)
AMREX_HOME = $(abspath ../../lib/amrex)
PELE_PHYSICS_HOME = $(abspath ../../lib/PelePhysics)

CEXE_headers += prob.h
include $(AMR_SOLVER)/src/Make.CNS
