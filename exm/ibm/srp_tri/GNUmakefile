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
PELE_PHYSICS_HOME = $(abspath ../../../lib/PelePhysics)

# Boost and CGAL locations (MacOS example)
BOOST_HOME=/opt/homebrew/opt/boost
CGAL_HOME=/opt/homebrew/opt/cgal
GMP_HOME=/opt/homebrew/opt/gmp
MPFR_HOME=/opt/homebrew/opt/mpfr

CEXE_headers += prob.h
include $(AMR_SOLVER)/src/Make.CNS
