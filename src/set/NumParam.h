#ifndef NUMPARAM_H_
#define NUMPARAM_H_

#include <AMReX_GpuContainers.H>

struct Num_Method_Param {

public:

bool dissipation = false;

#ifdef AMREX_USE_GPIBM
bool ibm = true;
#else
bool ibm = false;
#endif

int order = 2;                        // order of scheme
amrex::Real C2skew=0.1,C4skew=0.0016; // Skew symmetric default
amrex::Real C1=0.0,C2=0.0;            // dummy paremeters

};



#endif

