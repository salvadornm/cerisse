#ifndef CNS_PROB_PARM_H
#define CNS_PROB_PARM_H

#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>

struct ProbParm
{
  amrex::Real mach = 0.05;
  amrex::Real beta = 0.02; // vortex strength
  amrex::Real p0 = 1e6; // [erg cm^-3]
  amrex::Real T0 = 300.0; // [K]
  amrex::Real rho0;
  amrex::Real v0;
};

#endif
