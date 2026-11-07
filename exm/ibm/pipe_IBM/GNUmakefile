DEBUG     = FALSE
DIM       = 3
COMP      = gnu
CXXSTD    = c++20
PRECISION = DOUBLE

TINY_PROFILE = FALSE
CLIP_MINTEMP = FALSE

USE_MPI  = FALSE
USE_CUDA = FALSE
USE_OMP  = FALSE
USE_GPIBM = TRUE
USE_EB    = FALSE

USE_PELEPHYSICS = TRUE
EOS_MODEL := FUEGO
TRANSPORT_MODEL := SIMPLE
CHEMISTRY_MODEL := BurkeDryer
ENSURE_MASSFRACSUM = TRUE

AMREX_NO_PROBINIT = TRUE

AMR_SOLVER = $(abspath ../../../)
AMREX_HOME = $(abspath ../../../lib/amrex)
PELE_PHYSICS_HOME = $(abspath ../../../lib/PelePhysics)
SUNDIALS_HOME ?= $(abspath ../../../lib/sundials)

CEXE_headers += prob.h
include $(AMR_SOLVER)/src/Make.CNS
