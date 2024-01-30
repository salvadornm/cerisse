#ifndef CNS_PROB_PARM_H
#define CNS_PROB_PARM_H

#include <AMReX_GpuMemory.H>
#include <AMReX_REAL.H>
#include <math.h>

using namespace amrex::literals;

struct ProbParm
{
  amrex::Real Ma_oo=6.0;  // free-stream mach Number
  amrex::Real p_l = 5000; // free-stream pressure (at 10 km altitude)
  amrex::Real u_l;
  amrex::Real T_l = 223;  // free-stream temperture
  amrex::Real rho_l;
  amrex::Real rhoe_l;

};

#endif
