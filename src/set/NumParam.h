#ifndef NUMPARAM_H_
#define NUMPARAM_H_

#include <AMReX_GpuContainers.H>

using namespace amrex;

// default values 
struct defaultparm_t {

public:

static constexpr bool dissipation = false;

#ifdef AMREX_USE_GPIBM
static constexpr bool ibm = true;
#else
static constexpr bool ibm = false;
#endif

static constexpr int order = 4;                   // order of scheme
static constexpr Real C2skew=0.1,C4skew=0.0016;   // Skew symmetric default
static constexpr Real conductivity = 0.0262;      // conductivity (for constant value)
static constexpr Real viscosity   = 1.85e-5;      // viscosity    (for constant value)

static constexpr bool solve_viscterms_only = false;


};



#endif

