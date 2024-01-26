#ifndef CNS_PROB_PARM_H
#define CNS_PROB_PARM_H

#include <AMReX_GpuMemory.H>
#include <AMReX_REAL.H>

struct ProbParm
{
  amrex::Real T0 = 800.0; // [K]
  amrex::Real p0 = 2.0e5; // [Ba]

  amrex::Real L;
  amrex::Real rho0;
  amrex::Real v0;

  amrex::GpuArray<amrex::Real, NUM_SPECIES> Y = {0.0};
};

#endif
