#ifndef CNS_PROB_PARM_H
#define CNS_PROB_PARM_H

#include <AMReX_GpuMemory.H>
#include <AMReX_REAL.H>
#include <math.h>

struct ProbParm
{
  amrex::Real rho;
  amrex::Real u;
  amrex::Real ei;
  amrex::GpuArray<amrex::Real, NUM_SPECIES> Y = {0.0};

  amrex::Real rho_j;
  amrex::Real v_j;
  amrex::Real ei_j;
  amrex::GpuArray<amrex::Real, NUM_SPECIES> Y_jet = {0.0};

  bool spark = false;
};

#endif