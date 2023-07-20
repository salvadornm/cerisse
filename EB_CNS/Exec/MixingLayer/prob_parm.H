#ifndef CNS_PROB_PARM_H
#define CNS_PROB_PARM_H

#include <AMReX_GpuMemory.H>
#include <AMReX_REAL.H>

struct ProbParm
{
  amrex::Real T1 = 545;
  amrex::Real u1 = 66910;
  amrex::Real rho1;
  amrex::Real ei1;
  amrex::GpuArray<amrex::Real, NUM_SPECIES> massfrac1 = {0.0};

  amrex::Real T2 = 1475;
  amrex::Real u2 = 115160;
  amrex::Real rho2;
  amrex::Real ei2;
  amrex::GpuArray<amrex::Real, NUM_SPECIES> massfrac2 = {0.0};

  amrex::Real p = 942322.5;
  amrex::Real theta_w = 0.0198; // vorticity thickness
  amrex::Real uc;               // convective velocity

  bool record_statistics; // record statistics or not
};

#endif
