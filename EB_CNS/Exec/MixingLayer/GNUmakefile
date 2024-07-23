# AMReX
DIM = 2
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
Custom_Chemistry_Model := ../../Source/chem/OConaire-9-21

# CNS
NUM_FIELD = 0
NUM_AUX = 13 # for time averages

# GNU Make
include ./Make.package
include ../Make.CNS