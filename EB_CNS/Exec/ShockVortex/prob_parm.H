#ifndef CNS_PROB_PARM_H
#define CNS_PROB_PARM_H

#include <AMReX_GpuMemory.H>
#include <AMReX_REAL.H>
#include <math.h>

using namespace amrex::literals;

struct ProbParm
{
  // Upstream
  amrex::Real p_l = 1.0;
  amrex::Real u_l = 1.5 * sqrt(1.4);
  amrex::Real rho_l = 1.0;
  amrex::Real rhoe_l = p_l / (1.4 - 1.0);

  // Downstream (depends on M)
  amrex::Real p_r;
  amrex::Real u_r;
  amrex::Real rho_r;
};

#endif
