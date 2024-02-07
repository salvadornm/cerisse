#ifndef CNS_PROB_PARM_H
#define CNS_PROB_PARM_H

#include <AMReX_GpuMemory.H>
#include <AMReX_REAL.H>

using namespace amrex::literals;

struct ProbParm
{
  amrex::Real p_l = 71730;      // left pressure [0.1*Pa]
  amrex::Real u_l = 0.0;        // left velocity [cm/s]
  amrex::Real rho_l = 0.072e-3; // left density [g/cm^3]
  amrex::Real e_l;
  amrex::GpuArray<amrex::Real, NUM_SPECIES> massfrac_l = {0.0};

  amrex::Real p_r = 355940;       // right pressure [0.1*Pa]
  amrex::Real u_r = -48734;       // right velocity [cm/s]
  amrex::Real rho_r = 0.18075e-3; // right density [g/cm^3]
  amrex::Real e_r;
  amrex::GpuArray<amrex::Real, NUM_SPECIES> massfrac_r = {0.0};
};

#endif
