#ifndef CNS_PROB_PARM_H
#define CNS_PROB_PARM_H

#include <AMReX_GpuMemory.H>
#include <AMReX_REAL.H>
#include <math.h>

using namespace amrex::literals;

struct ProbParm
{

  amrex::Real p_l = 153546.24; // post-shock
  amrex::Real u_l;
  // amrex::Real v_l = 0.0;
  amrex::Real T_l = 169.167;
  amrex::Real rho_l;
  amrex::Real rhoe_l;

  amrex::Real p_r = 153546.24; // undisturbed
  amrex::Real u_r = 0.0;
  // amrex::Real v_r = 0.0;
  amrex::Real T_r = 153.0;
  amrex::Real rho_r;
  amrex::Real rhoe_r;
  // amrex::Real c_r;
};

#endif
