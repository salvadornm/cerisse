#ifndef CNS_CONSTANTS_H
#define CNS_CONSTANTS_H

#include <AMReX_GpuMemory.H>
#include <AMReX_REAL.H>

namespace CNSConstants {

constexpr amrex::Real smallr   = 1.0e-19;   // small rho
constexpr amrex::Real smallp   = 1.0e-10;   // small p
constexpr amrex::Real smallu   = 1.0e-12;   // small u
constexpr amrex::Real smally   = 1.0e-10;   // tol of |sum(Y)-1|
constexpr amrex::Real smalleps = 1.0e-8;    // tol of 


constexpr int level_mask_interior = 0;    // valid cells
constexpr int level_mask_covered = 1;     // ghost cells covered by valid cells of this level
constexpr int level_mask_notcovered = 2;  // ghost cells not covered
constexpr int level_mask_physbnd = 3;     // outside domain

constexpr amrex::Real min_react_temp  = 350.0;   // minimum reaction temperature
constexpr amrex::Real min_euler_temp  = 50.0;    // minimum temperature for Euler solvers  TODO: make it a user parameter?
constexpr amrex::Real min_euler_press = 1.0e-8;  // minimum pressure for Euler solvers

constexpr amrex::Real one_half  = 1.0/2.0;   // 1/2
constexpr amrex::Real one_third = 1.0/3.0;   // 1/3
constexpr amrex::Real two_third = 2.0/3.0;   // 2/3

 
}; // namespace CNSConstants

#endif
