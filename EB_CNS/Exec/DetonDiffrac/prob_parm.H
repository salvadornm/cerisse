#ifndef CNS_PROB_PARM_H
#define CNS_PROB_PARM_H

#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>

using namespace amrex::literals;

struct ProbParm
{
    amrex::Real p_1 = 1.5; // left
    amrex::Real u_1 = 0.0;
    amrex::Real v_1 = 0.0;
    amrex::Real rho_1 = 1.5;
    amrex::Real rhoe_1;
    amrex::GpuArray<amrex::Real, NUM_SPECIES> massfrac = {0.0};

    amrex::Real p_2 = 0.3; // right
    amrex::Real u_2 = 1.206;
    amrex::Real v_2 = 0.0;
    amrex::Real rho_2 = 0.5323;
    amrex::Real rhoe_2;
    amrex::GpuArray<amrex::Real, NUM_SPECIES> massfrac_2 = {0.0};
};

#endif
