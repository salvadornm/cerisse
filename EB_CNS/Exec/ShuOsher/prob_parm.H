#ifndef CNS_PROB_PARM_H
#define CNS_PROB_PARM_H

#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>

using namespace amrex::literals;

struct ProbParm
{
    // amrex::Real p_l = 10.3333;      // left pressure
    // amrex::Real u_l = 2.629369;      // left velocity
    // amrex::Real rho_l = 3.857143; // left density
    // amrex::Real rhoe_l;
    // amrex::GpuArray<amrex::Real, NUM_SPECIES> massfrac_l = {1.0};
    // amrex::Real sd_u_l = 2000.0;  // velocity std dev

    // amrex::Real p_r = 1.0;      // right pressure
    // amrex::Real u_r = 0.0;      // right velocity
    // amrex::Real rho_r = 1.0; // right density
    // amrex::Real rhoe_r;
    // amrex::GpuArray<amrex::Real, NUM_SPECIES> massfrac_r = {1.0};
};

#endif
