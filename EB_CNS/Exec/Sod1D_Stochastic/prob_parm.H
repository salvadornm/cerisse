#ifndef CNS_PROB_PARM_H
#define CNS_PROB_PARM_H

#include <AMReX_GpuMemory.H>
#include <AMReX_REAL.H>

using namespace amrex::literals;

// struct ProbParm
// {
//     amrex::Real p_l = 1e6;   // left pressure
//     amrex::Real u_l = 0.0;   // left velocity
//     amrex::Real rho_l = 0.729e-3; // left density
//     amrex::Real rhoe_l;
//     amrex::Real T_l;
//     amrex::GpuArray<amrex::Real, NUM_SPECIES> massfrac_l;
//     amrex::Real sd_u_l = 50; // velocity std dev

//     amrex::Real p_r = 5e5;     // right pressure
//     amrex::Real u_r = 0.0;     // right velocity
//     amrex::Real rho_r = 0.456e-3; // right density
//     amrex::Real rhoe_r;
//     amrex::Real T_r;
//     amrex::GpuArray<amrex::Real, NUM_SPECIES> massfrac_r;
// };

struct ProbParm
{
  amrex::Real p_l = 1e7;      // left pressure
  amrex::Real u_l = 0.0;      // left velocity
  amrex::Real rho_l = 9.6e-4; // left density
  amrex::Real rhoe_l;
  amrex::GpuArray<amrex::Real, NUM_SPECIES> massfrac_l;
  amrex::Real sd_u_l = 2000.0; // velocity std dev

  amrex::Real p_r = 1e6;      // right pressure
  amrex::Real u_r = 0.0;      // right velocity
  amrex::Real rho_r = 1.2e-4; // right density
  amrex::Real rhoe_r;
  amrex::GpuArray<amrex::Real, NUM_SPECIES> massfrac_r;
};

#endif
