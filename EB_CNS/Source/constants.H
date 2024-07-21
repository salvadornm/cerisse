#ifndef CNS_CONSTANTS_H
#define CNS_CONSTANTS_H

#include <AMReX_GpuMemory.H>
#include <AMReX_REAL.H>

namespace CNSConstants {
constexpr amrex::Real smallr = 1.0e-19; // small rho
// const amrex::Real maxr = 100.0; // biggest rho
constexpr amrex::Real smallp = 1.0e-10; // small p
constexpr amrex::Real smallu = 1.0e-12; // small u
constexpr amrex::Real smally = 1.0e-10; // tol of |sum(Y)-1|

constexpr int level_mask_interior = 0; // valid cells
constexpr int level_mask_covered = 1; // ghost cells covered by valid cells of this level
constexpr int level_mask_notcovered = 2; // ghost cells not covered
constexpr int level_mask_physbnd = 3;    // outside domain
}; // namespace CNSConstants

#endif
