# IBM test: Complex geometry (multiple polygons, BVH CPU 2D)
DEBUG     = FALSE
DIM       = 2
COMP      = gcc
CXXSTD    = c++20
PRECISION = DOUBLE
TINY_PROFILE = FALSE
CLIP_MINTEMP = FALSE

USE_MPI  = TRUE
USE_CUDA = FALSE
USE_OMP  = FALSE
USE_GPIBM = TRUE

AMREX_NO_PROBINIT = TRUE
CUDA_ARCH = 8.0

AMR_SOLVER = $(abspath ../../../../)
AMREX_HOME = $(AMR_SOLVER)/lib/amrex

# files in  TEST_DIR
TEST_DIR= $(AMR_SOLVER)/exm/ibm/ibm_tests/complex_geom
CEXE_headers += $(TEST_DIR)/prob.h
INCLUDE_LOCATIONS += $(TEST_DIR)

include $(AMR_SOLVER)/src/Make.CNS
