#ifndef NUMPARAM_H_
#define NUMPARAM_H_

#include <AMReX_GpuContainers.H>

using namespace amrex;

// default values 
struct defaultparm_t {

public:


static constexpr int order = 2;                   // order of scheme (used in many options)

static constexpr bool dissipation = false;        // use in some schemes to activate dissipation

#ifdef AMREX_USE_GPIBM
static constexpr bool ibm = true;
#else
static constexpr bool ibm = false;
#endif

//  Skew 
static constexpr Real C2skew=0.1,C4skew=0.0016;   // Skew symmetric default

// Transport properties
static constexpr Real conductivity = 0.0262;      // conductivity (for constant value)
static constexpr Real viscosity   = 1.85e-5;      // viscosity    (for constant value)

// ebm/walls options
static constexpr Real Twall = 300;              // wall temperature (for isothermal wall type)
static constexpr bool solve_diffwall = false;   // do not solve viscous fluxes at teh wall (for ebm)

// viscous options 
static constexpr bool use_LES = false;

// LES options
static constexpr Real Pr_o_Prsgs = 1.0;          // Pr/Prsgs
static constexpr Real Scsgs = 0.7;              // sgs Schmidt number
static constexpr Real Cs = 0.1;                 // Smagorinsky constant
static constexpr Real CI = 0.1;                 // Yoshizawa constant
static constexpr bool fixDelta = false;         // (it will use local lmesh size as filter width)
static constexpr Real Delta = 0.01;             // Filter width  L/20

};



#endif

