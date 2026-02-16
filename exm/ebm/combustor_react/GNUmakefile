DEBUG   = FALSE
DIM     = 3
COMP     = gnu
CXXSTD    = c++20
PRECISION = DOUBLE

TINY_PROFILE = FALSE

CLIP_MINTEMP = FALSE

USE_MPI = TRUE
USE_CUDA = FALSE
USE_OMP = FALSE

# Boundaries
USE_GPIBM=FALSE
USE_EB = TRUE


# PelePhysics
EOS_MODEL := FUEGO
TRANSPORT_MODEL := SIMPLE
CHEMISTRY_MODEL := LiDryer
USE_PELEPHYSICS = TRUE

# Cerisse Option ensure sumY=1
ENSURE_MASSFRACSUM = TRUE

# Cerisse Option
#CERIS_GPU = FALSE

AMREX_NO_PROBINIT=TRUE
#CUDA_ARCH=8.0

AMR_SOLVER = $(abspath ../../../)
AMREX_HOME = $(abspath ../../../lib/amrex)
PELE_PHYSICS_HOME = $(abspath ../../../lib/PelePhysics)
SUNDIALS_HOME = $(abspath ../../../lib/sundials/)
PELEPVERSION = 25
AMREXVERSION = 25


CEXE_headers += prob.h
CEXE_sources += eb_geometry.cpp
include $(AMR_SOLVER)/src/Make.CNS
