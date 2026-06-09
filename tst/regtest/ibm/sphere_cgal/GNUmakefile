DEBUG     = FALSE
DIM       = 3
COMP      = gcc
CXXSTD    = c++20
PRECISION = DOUBLE

TINY_PROFILE = FALSE

# cerisse options
CLIP_MINTEMP = FALSE

# parallelization
USE_MPI = TRUE
USE_CUDA = FALSE
USE_OMP = FALSE

# Boundaries
USE_GPIBM=TRUE
USE_CGAL  = TRUE

AMREX_NO_PROBINIT=TRUE
CUDA_ARCH=8.0

# Paths
AMR_SOLVER = $(abspath ../../../../)
AMREX_HOME = $(AMR_SOLVER)/lib/amrex


#  possible Boost and CGAL locations (searching homebrew first)
BOOST_PREFIX=/opt/homebrew/opt/boost
CGAL_PREFIX=/opt/homebrew/opt/cgal

ifeq ($(wildcard $(BOOST_PREFIX)/include/boost),)
  # Boost not found at Homebrew path — (default will be used)
else
  BOOST_HOME := $(BOOST_PREFIX)
  CGAL_HOME := $(CGAL_PREFIX)
  GMP_HOME=/opt/homebrew/opt/gmp
  MPFR_HOME=/opt/homebrew/opt/mpfr
endif

# IBM test: CGAL CPU 3D (single sphere)
# files in  TEST_DIR
TEST_DIR= $(AMR_SOLVER)/exm/ibm/ibm_tests/
CEXE_headers += $(TEST_DIR)/3d_cgal/prob.h
INCLUDE_LOCATIONS += $(TEST_DIR)
VPATH_LOCATIONS   +=  $(TEST_DIR)
#

include $(AMR_SOLVER)/src/Make.CNS
