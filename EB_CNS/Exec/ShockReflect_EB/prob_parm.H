#ifndef CNS_PROB_PARM_H
#define CNS_PROB_PARM_H

#include <AMReX_GpuMemory.H>
#include <AMReX_REAL.H>
#include <math.h>

using namespace amrex::literals;

struct ProbParm
{
  amrex::Real p_l = 116.5; // post-shock
  amrex::Real u_l = 8.25;
  amrex::Real v_l = 0.0;
  amrex::Real rho_l = 8.0;
  amrex::Real rhoe_l;
  // amrex::Real T_l;

  amrex::Real p_r = 1.0; // undisturbed
  amrex::Real u_r = 0.0;
  amrex::Real v_r = 0.0;
  amrex::Real rho_r = 1.4;
  amrex::Real rhoe_r;
  // amrex::Real T_r;
  amrex::Real c_r;
};

#endif
