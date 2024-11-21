TINY_PROFILE = FALSE

DEBUG   = FALSE
DIM     = 2
COMP     = gcc
CXXSTD    = c++20

USE_MPI = TRUE
USE_CUDA = FALSE
USE_OMP = FALSE
USE_GPIBM=FALSE


# PelePhysics
USE_PELEPHYSICS = FALSE
#EOS_MODEL := FUEGO
#TRANSPORT_MODEL := SIMPLE
#CHEMISTRY_MODEL := BurkeDryer

AMREX_NO_PROBINIT=TRUE
CUDA_ARCH=8.0

# Paths
AMR_SOLVER = $(abspath ../../)
AMREX_HOME = $(abspath ../../lib/amrex)
PELE_PHYSICS_HOME = $(abspath ../../lib/PelePhysics)

CEXE_headers += prob.h
include $(AMR_SOLVER)/src/Make.CNS
