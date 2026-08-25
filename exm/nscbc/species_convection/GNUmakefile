DIM = 2
COMP = gnu
CXXSTD = c++20
PRECISION = DOUBLE

TINY_PROFILE = FALSE
DEBUG = FALSE
USE_MPI = TRUE
USE_OMP = FALSE

USE_PELEPHYSICS = TRUE
EOS_MODEL := FUEGO
TRANSPORT_MODEL := SIMPLE
CHEMISTRY_MODEL := air

AMREX_NO_PROBINIT = TRUE
USE_GPIBM = FALSE
USE_EB = FALSE

AMR_SOLVER = $(abspath ../../../)
AMREX_HOME = $(abspath ../../../lib/amrex)
PELE_PHYSICS_HOME = $(abspath ../../../lib/PelePhysics)
SUNDIALS_HOME = $(abspath ../../../lib/sundials)

PELEPVERSION = 25
AMREXVERSION = 25

include ./Make.package
include $(AMR_SOLVER)/src/Make.CNS
