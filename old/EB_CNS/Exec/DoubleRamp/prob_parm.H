#ifndef CNS_PROB_PARM_H
#define CNS_PROB_PARM_H

#include <AMReX_GpuMemory.H>
#include <AMReX_REAL.H>
#include <math.h>

using namespace amrex::literals;

struct ProbParm
{
  amrex::Real M = 3.0;
  amrex::Real p = 204170.0;
  amrex::Real T = 103.0;
  amrex::Real rho;
  amrex::Real u;
  amrex::Real rhoe;
  amrex::GpuArray<amrex::Real, NUM_SPECIES> massfrac = {0.0};
};

#endif
