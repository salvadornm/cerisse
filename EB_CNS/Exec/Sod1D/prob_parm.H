#ifndef CNS_PROB_PARM_H
#define CNS_PROB_PARM_H

#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>

using namespace amrex::literals;

struct ProbParm
{
    amrex::Real p_l = 1.0;   // left pressure
    amrex::Real u_l = 0.0;   // left velocity
    amrex::Real rho_l = 1.0; // left density
    amrex::Real rhoe_l;
    amrex::Real T_l;

    amrex::Real p_r = 0.1;     // right pressure
    amrex::Real u_r = 0.0;     // right velocity
    amrex::Real rho_r = 0.125; // right density
    amrex::Real rhoe_r;
    amrex::Real T_r;
};

#endif
