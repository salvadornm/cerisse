DEBUG   = FALSE
DIM     = 3
COMP     = gcc
CXXSTD    = c++23
PRECISION = DOUBLE

TINY_PROFILE = FALSE

# cerisse options
CLIP_MINTEMP = FALSE

# parallelization
USE_MPI = TRUE
USE_CUDA = FALSE
USE_OMP = FALSE

# PelePhysics
USE_PELEPHYSICS = FALSE

# Boundaries
USE_GPIBM=TRUE

AMREX_NO_PROBINIT=TRUE
CUDA_ARCH=8.0

# Paths
AMR_SOLVER = $(abspath ../../../)
AMREX_HOME = $(abspath ../../../lib/amrex)

# Boost and CGAL locations (Local installation)
# GMP and MPFR (default)
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

CEXE_headers += prob.h
include $(AMR_SOLVER)/src/Make.CNS
