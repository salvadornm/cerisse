TINY_PROFILE 	= FALSE
DEBUG   	    = FALSE
PRECISION     = DOUBLE # DOUBLE or FLOAT
DIM     	    = 3
COMP 		      = gcc
CXXSTD        = c++17 # c++17 default
USE_CUDA 	= FALSE
USE_MPI 	= FALSE
USE_OMP 	= FALSE
USE_GPIBM	= TRUE
AMREX_NO_PROBINIT =TRUE

include ./Make.package
include $(AMREX_HOME)/../../src/Make.CNS

INCLUDE_LOCATIONS += $(AMREX_HOME)/../install/boost/include
LIBRARY_LOCATIONS += $(AMREX_HOME)/../install/boost/lib

INCLUDE_LOCATIONS += $(AMREX_HOME)/../install/cgal/include
LIBRARY_LOCATIONS += $(AMREX_HOME)/../install/cgal/lib/cmake/CGAL
LIBRARIES += -lgmp