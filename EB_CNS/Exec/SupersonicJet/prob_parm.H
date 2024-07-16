#ifndef CNS_PROB_PARM_H
#define CNS_PROB_PARM_H

#include <AMReX_GpuMemory.H>
#include <AMReX_REAL.H>

struct ProbParm
{
  // Vitiated air inlet
  amrex::Real T1 = 1350; //1250;
  amrex::Real u1 = 142000;
  amrex::Real p1 = 1070000;
  amrex::Real rho1;
  amrex::Real ei1;
  amrex::GpuArray<amrex::Real, NUM_SPECIES> Y1 = {0.0};

  // Fuel inlet
  amrex::Real T2 = 545;
  amrex::Real u2 = 178000;
  amrex::Real p2 = 1120000;
  amrex::Real rho2;
  amrex::Real ei2;
  amrex::GpuArray<amrex::Real, NUM_SPECIES> Y2 = {0.0};

  // Atmospheric
  amrex::Real T_a = 296;
  amrex::Real p_a = 1010000;
  amrex::Real rho_a;
  amrex::Real ei_a;
  amrex::GpuArray<amrex::Real, NUM_SPECIES> Y_a = {0.0};

  bool record_statistics = false;
  bool clean_aux_on_restart = false;
  bool inflow_turbulence = false;
};

#endif
