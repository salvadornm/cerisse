# AMReX
DIM = 3
COMP = gnu
PRECISION = DOUBLE

# Profiling
PROFILE = FALSE
TINY_PROFILE = FALSE
COMM_PROFILE = FALSE
TRACE_PROFILE = FALSE
MEM_PROFILE = FALSE
USE_GPROF = FALSE

# Performance
USE_MPI = TRUE
USE_OMP = FALSE
USE_CUDA = FALSE
USE_HIP = FALSE
USE_SYCL = FALSE

# Debugging
DEBUG = FALSE
FSANITIZER = FALSE
THREAD_SANITIZER = FALSE

# PelePhysics
Eos_Model := Fuego
Transport_Model := Simple
#Chemistry_Model := LiDryer
Custom_Chemistry_Model := FakeReaction2

# CNS
USE_FULL_PROB_POST_TIMESTEP = TRUE
USE_PROB_PARM_HOST = TRUE
NUM_FIELD = 8

# GNU Make
CEXE_headers += prob.H prob_parm.H
CEXE_sources += prob.cpp
include ../../Source/Make.CNS
