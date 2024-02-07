#ifndef CNS_PROB_PARM_H
#define CNS_PROB_PARM_H

#include <AMReX_GpuMemory.H>
#include <AMReX_REAL.H>

using namespace amrex::literals;

struct ProbParm
{
  amrex::Real p0 = 1.01325e6;
  amrex::Real x_offset;
  amrex::Real u_offset;
};

#endif
