#ifndef CNS_CONSTANTS_H
#define CNS_CONSTANTS_H

using namespace amrex;

namespace CNSConstants {

    constexpr Real smallr   = 1.0e-19;   // small rho
    constexpr Real smallp   = 1.0e-10;   // small p
    constexpr Real smallu   = 1.0e-12;   // small u
    constexpr Real smally   = 1.0e-10;   // tol of |sum(Y)-1|
    constexpr Real smalleps = 1.0e-8;    // tol of 

    constexpr int level_mask_interior = 0;    // valid cells
    constexpr int level_mask_covered = 1;     // ghost cells covered by valid cells of this level
    constexpr int level_mask_notcovered = 2;  // ghost cells not covered
    constexpr int level_mask_physbnd = 3;     // outside domain

    constexpr Real min_react_temp  = 350.0;   // minimum reaction temperature
    constexpr Real min_euler_temp  = 10.0;    // minimum temperature for Euler solvers
    constexpr Real min_euler_press = 1.0e-8;  // minimum pressure for Euler solvers
    
    // Tiny functions to avoid reference  device vars (better for GPU)
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static constexpr Real min_press() noexcept { return Real(min_euler_press); }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static constexpr Real min_temp() noexcept { return Real(min_euler_temp); }
    //

    constexpr Real one_half  = 1.0/2.0;   // 1/2
    constexpr Real one_third = 1.0/3.0;   // 1/3
    constexpr Real two_third = 2.0/3.0;   // 2/3
 
}; // namespace CNSConstants

#endif
