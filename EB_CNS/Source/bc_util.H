#ifndef CNS_BC_UTIL_H_
#define CNS_BC_UTIL_H_

#include <AMReX_REAL.H>

#include "index_macros.H"

/**
 * @brief Helper function to generate ghost cell states for adiabatic wall.
 *
 * @param[out] s_ext output ghost cell states.
 * @param s_int    input near-wall states.
 * @param dir      direction, 0:x, 1:y, 2:z.
 * @param no_slip  no-slip or not
 */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
bc_adiabatic_wall(amrex::Real s_ext[NVAR], const amrex::Real s_int[NVAR], int dir,
                  bool no_slip = true)
{
  s_ext[URHO] = s_int[URHO];
  s_ext[UMX] = ((dir == 0 || no_slip) ? -1 : 1) * s_int[UMX];
  s_ext[UMY] = ((dir == 1 || no_slip) ? -1 : 1) * s_int[UMY];
  s_ext[UMZ] = ((dir == 2 || no_slip) ? -1 : 1) * s_int[UMZ];
  s_ext[UEDEN] = s_int[UEDEN];
  for (int n = 0; n < NUM_SPECIES; ++n) s_ext[UFS + n] = s_int[UFS + n];
}

/**
 * @brief Helper function to generate ghost cell states for isothermal wall.
 *
 * @param[out] s_ext output ghost cell states.
 * @param s_int    input near-wall states.
 * @param Tw       wall temperature.
 * @param dir      direction, 0:x, 1:y, 2:z.
 * @param no_slip  no-slip or not (default true)
 */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
bc_isothermal_wall(amrex::Real s_ext[NVAR], const amrex::Real s_int[NVAR],
                   amrex::Real Tw, int dir, bool no_slip = true)
{
  auto eos = pele::physics::PhysicsType::eos();

  // Fluid states
  amrex::Real rhoinv = 1.0 / s_int[URHO];
  amrex::Real y[NUM_SPECIES], T, p;
  for (int n = 0; n < NUM_SPECIES; ++n) y[n] = s_int[UFS + n] * rhoinv;
  amrex::Real ke = 0.5 * rhoinv * rhoinv *
                   (AMREX_D_TERM(s_int[UMX] * s_int[UMX], +s_int[UMY] * s_int[UMY],
                                 +s_int[UMZ] * s_int[UMZ]));
  amrex::Real ei = rhoinv * s_int[UEDEN] - ke;
  eos.REY2T(s_int[URHO], ei, y, T);
  eos.RTY2P(s_int[URHO], T, y, p);

  // Wall states (extrapolate p and y, calculate rhow, ew based on Tw)
  amrex::Real rhow, ew;
  amrex::Real Tg =
    amrex::max(1e-5, 2 * Tw - T); // temp of ghost cell, prevent negative temp
  eos.PYT2RE(p, y, Tg, rhow, ew);
  ew += ke;

  s_ext[URHO] = rhow;
  s_ext[UMX] = ((dir == 0 || no_slip) ? -1 : 1) * rhow * rhoinv * s_int[UMX];
  s_ext[UMY] = ((dir == 1 || no_slip) ? -1 : 1) * rhow * rhoinv * s_int[UMY];
  s_ext[UMZ] = ((dir == 2 || no_slip) ? -1 : 1) * rhow * rhoinv * s_int[UMZ];
  s_ext[UEDEN] = rhow * ew;
  for (int n = 0; n < NUM_SPECIES; ++n) s_ext[UFS + n] = rhow * y[n];
}

// TODO: Riemann invariant based in/outlet, NSCBC

// /**
//  * @brief Helper function to generate ghost cell states for constant
//  momentum-flow-rate inflow.
//  *
//  * @param[out] s_ext output ghost cell states.
//  * @param s_int    input near-wall states.
//  * @param rho      inflow density (only used when inflow becomes supersonic).
//  * @param T0       inflow total temperature.
//  * @param mx,my,mz inflow momentum.
//  * @param Y_in     inflow species mass fraction.
//  * @param dir      direction, 0:x, 1:y, 2:z.
//  */
// AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
// bc_const_mom_inlet (amrex::Real s_ext[NVAR],
//                     const amrex::Real s_int[NVAR],
//                     amrex::Real rho,
//                     amrex::Real T0,
//                     amrex::Real mx, amrex::Real my, amrex::Real mz,
//                     amrex::Real Y[NUM_SPECIES],
//                     int dir)
// {
//   // Extrapolate rho, Y, use user-specified momentum
//   amrex::Real rho_in = s_int[URHO];

//   auto eos = pele::physics::PhysicsType::eos();
//   amrex::Real e_in;
//   eos.RTY2Ei(rho_in, Y, T0, e_in); // calculate total energy directly

//   // Check if inflow is subsonic and temperature is > 0
//   amrex::Real T, cs;
//   amrex::Real vsqr = (mx*mx + my*my + mz*mz)/(rho_in*rho_in);
//   amrex::Real ei = e_in - 0.5*vsqr;
//   eos.REY2T(rho_in, ei, Y, T);
//   AMREX_ALWAYS_ASSERT(T > 0);

//   eos.RTY2Cs(rho_in, T, Y, cs);
//   amrex::Real M = sqrt(vsqr) / cs;
//   if (M > 1.0) {
//     rho_in = rho;
//     eos.RTY2Ei(rho_in, Y, T0, e_in);

//     vsqr = (mx*mx + my*my + mz*mz)/(rho_in*rho_in);
//     ei = e_in - 0.5*vsqr;
//     eos.REY2T(rho_in, ei, Y, T);
//     AMREX_ALWAYS_ASSERT(T > 0);
//   }

//   s_ext[URHO] = rho_in;
//   s_ext[UMX] = mx;
//   s_ext[UMY] = my;
//   s_ext[UMZ] = mz;
//   s_ext[UEDEN] = rho_in * e_in;
//   for (int n = 0; n < NUM_SPECIES; ++n)
//     s_ext[UFS+n] = rho_in * Y[n];
// }

#endif