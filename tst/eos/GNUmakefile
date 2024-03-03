DIM := 3
USE_EB := FALSE
LAZY := TRUE
BL_NO_FORT = TRUE

AMR_SOLVER = $(abspath ../../)
AMREX_HOME = $(abspath ../../lib/amrex)
PELE_PHYSICS_HOME = $(abspath ../../lib/PelePhysics)

CEXE_sources += main.cpp
CEXE_headers += catch.hpp

###########################################################

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

# include $(AMR_SOLVER)/src/Make.package
# INCLUDE_LOCATIONS += $(AMR_SOLVER)/src
# VPATH_LOCATIONS   += $(AMR_SOLVER)/src

# include $(AMR_SOLVER)/src/set/Make.package
CEXE_headers += Index.h
INCLUDE_LOCATIONS += $(AMR_SOLVER)/src/set/
# VPATH_LOCATIONS   += $(AMR_SOLVER)/src/set/

# include $(AMR_SOLVER)/src/tim/Make.package
# INCLUDE_LOCATIONS += $(AMR_SOLVER)/src/tim/
# VPATH_LOCATIONS   += $(AMR_SOLVER)/src/tim/

# include $(AMR_SOLVER)/src/rhs/Make.package
# INCLUDE_LOCATIONS += $(AMR_SOLVER)/src/rhs/
# VPATH_LOCATIONS   += $(AMR_SOLVER)/src/rhs/

include $(AMR_SOLVER)/src/cls/Make.package
INCLUDE_LOCATIONS += $(AMR_SOLVER)/src/cls/
VPATH_LOCATIONS   += $(AMR_SOLVER)/src/cls/

CHEM_HOME = Support/Mechanism/Models/Null
  
Pdirs := Source Eos Transport $(CHEM_HOME)
Ppack += $(foreach dir, $(Pdirs), $(PELE_PHYSICS_HOME)/$(dir)/Make.package)
Plocs += $(foreach dir, $(Pdirs), $(PELE_PHYSICS_HOME)/$(dir))

include $(Ppack)
INCLUDE_LOCATIONS += $(Plocs)
VPATH_LOCATIONS   += $(Plocs)

# include $(PELE_PHYSICS_HOME)/ThirdParty/Make.ThirdParty
# INCLUDE_LOCATIONS += $(PELE_PHYSICS_HOME)/ThirdParty
# VPATH_LOCATIONS   += $(PELE_PHYSICS_HOME)/ThirdParty  

DEFINES += -DUSE_PELEPHYSICS
DEFINES += -DUSE_GAMMALAW_EOS
DEFINES += -DUSE_CONSTANT_TRANSPORT

include $(AMREX_HOME)/Src/Base/Make.package
# include $(AMREX_HOME)/Src/Boundary/Make.package
# include $(AMREX_HOME)/Src/AmrCore/Make.package
# include $(AMREX_HOME)/Src/Amr/Make.package
include $(AMREX_HOME)/Tools/GNUMake/Make.rules

# stop on first compile error
warning_flags +=-Wfatal-errors

# CNS uses a coarse grained OMP approach
DEFINES += -DAMREX_CRSEGRNDOMP

#
# For debug, see amrex/Tools/GNUMake/README.md
#