#ifndef CNS_PROB_PARM_H
#define CNS_PROB_PARM_H

#include <AMReX_GpuMemory.H>
#include <AMReX_REAL.H>
#include <math.h>

using namespace amrex::literals;

struct ProbParm
{
  amrex::Real p = 1089877.72748;
  amrex::Real T = 255.56;
  amrex::Real rho;
  amrex::Real u;
  amrex::Real v;
  amrex::Real rhoe;
};

#endif
