# CPU-only numerical regression tests; no CERISSE/PelePhysics executable needed.
AMREX_HOME := $(abspath ../../lib/amrex)
DEBUG = TRUE
DIM = 3
COMP = gnu
CXXSTD = c++20
USE_MPI = FALSE
USE_OMP = FALSE
USE_CUDA = FALSE
BL_NO_FORT = TRUE
TINY_PROFILE = FALSE
KEEP_TEST_SPECIES ?= 3
KEEP_TEST_BOUNDARY ?= NONE
include $(AMREX_HOME)/Tools/GNUMake/Make.defs
CEXE_sources += main.cpp
CEXE_headers += Keep.h CentralKEEP_SPEC.h
INCLUDE_LOCATIONS += $(abspath ../../src/rhs)
VPATH_LOCATIONS += $(abspath ../../src/rhs)
DEFINES += -DKEEP_TEST_SPECIES=$(KEEP_TEST_SPECIES)
ifeq ($(KEEP_TEST_BOUNDARY),EB)
DEFINES += -DCNS_USE_EB=1
else ifeq ($(KEEP_TEST_BOUNDARY),IBM)
DEFINES += -DAMREX_USE_GPIBM=1
endif
include $(AMREX_HOME)/Src/Base/Make.package
include $(AMREX_HOME)/Tools/GNUMake/Make.rules
