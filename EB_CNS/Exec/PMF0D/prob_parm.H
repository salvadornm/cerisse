#ifndef CNS_PROB_PARM_H
#define CNS_PROB_PARM_H

#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>

using namespace amrex::literals;

struct ProbParm
{
    amrex::Real rho;
    amrex::Real eint;
    amrex::GpuArray<amrex::Real, NUM_SPECIES> massfrac = {0.0};
};

#endif
