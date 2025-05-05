# AMReX
DIM = 1
COMP = gcc
CXXSTD    = c++20
PRECISION = DOUBLE

# Profiling
TINY_PROFILE = FALSE

# Performance
USE_MPI = FALSE
USE_OMP = FALSE

# Debugging
DEBUG = FALSE

# PelePhysics
EOS_MODEL := FUEGO
TRANSPORT_MODEL := SIMPLE
CHEMISTRY_MODEL := LiDryer

# CNS
AMREX_NO_PROBINIT = TRUE
USE_GPIBM = FALSE
USE_EB = FALSE
USE_PELEPHYSICS = TRUE

# GNU Make
AMR_SOLVER = $(abspath ../../)
AMREX_HOME = $(abspath ../../lib/amrex)
PELE_PHYSICS_HOME = $(abspath ../../lib/PelePhysics)

CEXE_headers += prob.h
include $(AMR_SOLVER)/src/Make.CNS
