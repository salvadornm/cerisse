#ifndef CNS_PROB_PARM_H
#define CNS_PROB_PARM_H

#include <AMReX_GpuMemory.H>
#include <AMReX_REAL.H>

#include "mechanism.H"

using namespace amrex::literals;

struct ProbParm
{
  amrex::Real inflow_T = 285.5;    // [K]
  amrex::Real inflow_p = 1013250;    // [cgs] 
  amrex::Real inflow_mach = 0.022;
  amrex::Real interior_T = 285.5;  // [K] 
  amrex::Real interior_p = 1013250;  
  amrex::GpuArray<amrex::Real, NUM_SPECIES> massfrac = {1.0};

  amrex::Real* inflow_state = nullptr;

  ProbParm();
  ~ProbParm();
};

#endif
