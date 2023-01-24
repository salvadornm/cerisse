#ifndef CNS_PROB_PARM_H
#define CNS_PROB_PARM_H

#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>

using namespace amrex::literals;

struct ProbParm
{
    amrex::Real rho_l;
    amrex::Real u_l;
    amrex::Real rhoe_l;
    amrex::GpuArray<amrex::Real, NUM_SPECIES> massfrac_l = {0.0};

    amrex::Real rho_r;
    amrex::Real u_r;
    amrex::Real rhoe_r;
    amrex::GpuArray<amrex::Real, NUM_SPECIES> massfrac_r = {0.0};
};

#endif
