#ifndef CNS_PROB_PARM_H
#define CNS_PROB_PARM_H

#include <AMReX_GpuMemory.H>
#include <AMReX_REAL.H>
#include <math.h>

using namespace amrex::literals;

struct ProbParm
{
  amrex::Real p = 1530.0;
  amrex::Real T = 57.8;
  amrex::Real rho;
  amrex::Real u;
  amrex::Real rhoe;

  amrex::Real Tw = 300.0;
};

#endif
