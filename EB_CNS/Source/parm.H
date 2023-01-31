#ifndef CNS_PARM_H
#define CNS_PARM_H

#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>

struct Parm
{
    const amrex::Real smallr = 1.0e-19; // small rho
    const amrex::Real smallp = 1.0e-10; // small p
    const amrex::Real smallu = 1.0e-12; // small u
    const amrex::Real smally = 1.0e-10; // tol of |sum(Y)-1|

    static constexpr int level_mask_interior   = 0; // valid cells
    static constexpr int level_mask_covered    = 1; // ghost cells covered by valid cells of this level
    static constexpr int level_mask_notcovered = 2; // ghost cells not covered
    static constexpr int level_mask_physbnd    = 3; // outside domain
};

#endif
