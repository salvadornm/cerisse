#ifndef CNS_BC_UTIL_H_
#define CNS_BC_UTIL_H_

#include <AMReX_REAL.H>

#include "index_macros.H"
#include "nscbc.H"

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
  s_ext[UTEMP] = s_int[UTEMP];
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
  amrex::Real y[NUM_SPECIES], T = Tw, p;
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
  s_ext[UTEMP] = Tw;
}

/**
 * @brief Helper function to generate states downstream a normal shock. 
 *        Assume constant gamma equals to that upstream.
 *
 * @param[out] s_ext output downstream states.
 * @param s_int      input upstream states.
 * @param Ms         shock Mach number (if Ms = M upstream, the shock will be stationary).
 * @param dir        direction, 0:x, 1:y, 2:z.
 */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
normal_shock_relations(amrex::Real s_ext[NVAR], const amrex::Real s_int[NVAR],
                       amrex::Real Ms, int dir = 0)
{
  auto eos = pele::physics::PhysicsType::eos();

  // Upstream states
  amrex::Real rhoinv = 1.0 / s_int[URHO];
  amrex::Real y[NUM_SPECIES], T = 0.0, p, g;
  for (int n = 0; n < NUM_SPECIES; ++n) y[n] = s_int[UFS + n] * rhoinv;
  amrex::Real ke = 0.5 * rhoinv * rhoinv *
                   (AMREX_D_TERM(s_int[UMX] * s_int[UMX], +s_int[UMY] * s_int[UMY],
                                 +s_int[UMZ] * s_int[UMZ]));
  amrex::Real ei = rhoinv * s_int[UEDEN] - ke;
  eos.RYET2P(s_int[URHO], y, ei, T, p);
  eos.RTY2G(s_int[URHO], T, y, g);

  // Shocked states
  s_ext[URHO] = s_int[URHO] * (g + 1.0) * Ms * Ms / (2.0 + (g - 1.0) * Ms * Ms);
  amrex::Real u_r = s_int[UMX + dir] / s_int[URHO] * (2.0 + (g - 1.0) * Ms * Ms) /
                    ((g + 1.0) * Ms * Ms);
  amrex::Real p_r = p * (1.0 + (2.0 * g / (g + 1.0)) * (Ms * Ms - 1.0));

  // Convert back to conservative states
  s_ext[UMX] = (dir == 0) ? s_ext[URHO] * u_r : s_int[UMX];
  s_ext[UMY] = (dir == 1) ? s_ext[URHO] * u_r : s_int[UMY];
  s_ext[UMZ] = (dir == 2) ? s_ext[URHO] * u_r : s_int[UMZ];
  amrex::Real ei_r;
  eos.RYP2E(s_ext[URHO], y, p_r, ei_r);
  s_ext[UEDEN] = s_ext[URHO] * ei_r +
                 0.5 / s_ext[URHO] *
                   (AMREX_D_TERM(s_ext[UMX] * s_ext[UMX], +s_ext[UMY] * s_ext[UMY],
                                 +s_ext[UMZ] * s_ext[UMZ]));
  for (int n = 0; n < NUM_SPECIES; ++n) s_ext[UFS + n] = s_ext[URHO] * y[n];
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
bc_foextrap_outlet(const int dir, const int lo_hi, amrex::Real s_ext[NVAR],
                   const amrex::Real s_int[NVAR])
{
  // FOExtrap
  s_ext[URHO] = s_int[URHO];
  AMREX_D_TERM(s_ext[UMX] = s_int[UMX];, s_ext[UMY] = s_int[UMY];
               , s_ext[UMZ] = s_int[UMZ];)
  s_ext[UEDEN] = s_int[UEDEN];
  for (int n = 0; n < NUM_SPECIES; n++) s_ext[UFS + n] = s_int[UFS + n];

  // Outflow: prohibit backflow
  if (s_ext[UMX + dir] * lo_hi > 0.0) {
    s_ext[UMX + dir] = 0.0;
    s_ext[UEDEN] -= 0.5 * s_int[UMX + dir] * s_int[UMX + dir] / s_int[URHO];
  }
}

//
// AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
// bc_pressure_outlet(const amrex::Real x[AMREX_SPACEDIM], amrex::Real s_ext[NVAR], const amrex::Real s_int[NVAR],
//                    const amrex::Real p_ext, amrex::GeometryData const& geomdata)
// {
//   // Initialize all quantities with FOExtrap values
//   amrex::Real rho = s_int[URHO];
//   AMREX_D_TERM(amrex::Real u = s_int[UMX] / rho;,
//                amrex::Real v = s_int[UMY] / rho;,
//                amrex::Real w = s_int[UMZ] / rho;)
//   amrex::Real eint = s_int[UEDEN] / rho - 0.5 * (AMREX_D_TERM(u * u, + v * v, + w * w));
//   amrex::Real massfrac[NUM_SPECIES];
//   for (int n = 0; n < NUM_SPECIES; n++) { massfrac[n] = s_int[UFS + n] / rho; }
//   amrex::Real T = 0.0;
//   auto eos = pele::physics::PhysicsType::eos();
//   eos.REY2T(rho, eint, massfrac, T);

//   // Following Blazek p 279, eq. 8.23 - reimplemented from EB-C10
//   // NOTE: it is assumed that mass fractions are constant across the boundary
//   // Interior state (point d)
//   const amrex::Real* prob_hi = geomdata.ProbHi();
//   const amrex::Real* dx = geomdata.CellSize();
//   const amrex::Real xd = prob_hi[0] - 0.5 * dx[0];
//   const amrex::Real rho_inter = rho;
//   AMREX_D_TERM(const amrex::Real u_inter = u;,
//                const amrex::Real v_inter = v;,
//                const amrex::Real w_inter = w;)
//   const amrex::Real T_inter = T;
//   amrex::Real p_inter = 0.0, cs_inter = 0.0;
//   eos.RTY2P(rho_inter, T_inter, massfrac, p_inter);
//   eos.RTY2Cs(rho_inter, T_inter, massfrac, cs_inter);

//   // Boundary state (point b)
//   const amrex::Real xb = prob_hi[0];
//   const amrex::Real pb = p_ext;
//   const amrex::Real rhob = rho_inter + (pb - p_inter) / (cs_inter * cs_inter);
//   AMREX_D_TERM(
//     const amrex::Real ub = u_inter + (p_inter - pb) / (rho_inter * cs_inter);,
//     const amrex::Real vb = v_inter;,
//     const amrex::Real wb = w_inter;)

//   // Ghost state (point a). Linear extrapolation from d and b
//   rho = (rhob - rho_inter) / (xb - xd) * (x[0] - xd) + rho_inter;
//   amrex::Real p = (pb - p_inter) / (xb - xd) * (x[0] - xd) + p_inter;
//   p = amrex::max(p, 0.0001 * p_inter);
//   AMREX_ALWAYS_ASSERT(std::isgreater(p, 0.0));
//   AMREX_ALWAYS_ASSERT(std::isgreater(rho, 0.0));

//   eos.RYP2E(rho, massfrac, p, eint);
//   AMREX_D_TERM(
//     u = (ub - u_inter) / (xb - xd) * (x[0] - xd) + u_inter;,
//     v = (vb - v_inter) / (xb - xd) * (x[0] - xd) + v_inter;,
//     w = (wb - w_inter) / (xb - xd) * (x[0] - xd) + w_inter;)

//   // Write to s_ext
//   s_ext[URHO] = rho;
//   AMREX_D_TERM(s_ext[UMX] = rho * u;,
//                s_ext[UMY] = rho * v;,
//                s_ext[UMZ] = rho * w;)
//   s_ext[UEDEN] = rho * (eint + 0.5 * (AMREX_D_TERM(u * u, + v * v, + w * w)));
//   for (int n = 0; n < NUM_SPECIES; n++) { s_ext[UFS + n] = rho * massfrac[n]; }
// }

// AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
// bl_profile(const amrex::Real y, const amrex::Real z, ProbParm const& pp, 
//            amrex::Real rho, amrex::Real u, amrex::Real ei) {
//   // Make boundary layer
//   amrex::Real d_wall = amrex::min(y, 2.54 - y); // dist from upper and lower walls
//   d_wall = amrex::min(d_wall, 1.905-amrex::Math::abs(z)); // dist from upper and lower walls
//   d_wall /= (0.1245 * 2.0); // d_wall/D
//   amrex::Real u_profile = d_wall > 0.0 ? tanh(5.0 * d_wall) : 0.0; // u/u_inf
//   u = pp.u * u_profile;

//   amrex::Real T_inf, p_inf;
//   amrex::Real Y[NUM_SPECIES];
//   for (int n = 0; n < NUM_SPECIES; ++n) { Y[n] = pp.Y[n]; }
//   auto eos = pele::physics::PhysicsType::eos();
//   ei = pp.ei;
//   eos.RYET2P(pp.rho, Y, ei, T_inf, p_inf); // get T, p_inf

//   amrex::Real T_profile = d_wall > 0.0 ? 1.0 + (0.5 * pp.u * pp.u) / pp.cv_Tinf * (1.0 - u_profile * u_profile) : 1.0; // T/T_inf

//   amrex::Real cv, cv_inf = pp.cv_Tinf / T_inf;
//   eos.TY2Cv(T_inf * T_profile, Y, cv);
//   T_profile *= cv_inf / cv;

//   rho = pp.rho / T_profile;
//   eos.RYP2E(rho, Y, p_inf, ei); // rho, p_inf -> ei
// }

// /**
//  * @brief Helper function to generate ghost cell states for subsonic pressure outlet.
//  *        It returns to FOExtrap when local Mach number >= 1.
//  *
//  * @param[out] s_ext output ghost cell states.
//  * @param s_int    input near-wall states.
//  * @param p_ext    far field pressure.
//  * @param dir      direction, 0:x, 1:y, 2:z.
//  */
// AMREX_GPU_DEVICE AMREX_FORCE_INLINE void bc_pressure_outlet(
//   amrex::Real s_ext[NVAR], const amrex::Real s_int[NVAR], amrex::Real p_ext, int dir)
// {
//   auto eos = pele::physics::PhysicsType::eos();

//   // Fluid states
//   amrex::Real rhoinv = 1.0 / s_int[URHO];
//   amrex::Real y[NUM_SPECIES], T = 0.0, cs;
//   for (int n = 0; n < NUM_SPECIES; ++n) y[n] = s_int[UFS + n] * rhoinv;
//   amrex::Real ke = 0.5 * rhoinv * rhoinv *
//                    (AMREX_D_TERM(s_int[UMX] * s_int[UMX], +s_int[UMY] * s_int[UMY],
//                                  +s_int[UMZ] * s_int[UMZ]));
//   amrex::Real ei = rhoinv * s_int[UEDEN] - ke;
//   eos.REY2T(s_int[URHO], ei, y, T);
//   eos.RTY2Cs(s_int[URHO], T, y, cs);
//   amrex::Real M = std::sqrt(2 * ke) / cs; // lcoal Mach number

//   // Ghost states
//   if (M >= 1) {
//     for (int n = 0; n < NVAR; ++n) s_ext[n] = s_int[n]; // FOExtrap
//   } else {
//     amrex::Real gamma, p_int, cs_ext;
//     eos.RTY2G(s_int[URHO], T, y, gamma);
//     eos.RTY2P(s_int[URHO], T, y, p_int);
//     amrex::Real rhog = s_int[URHO] * std::pow(p_ext / p_int, 1.0 / gamma);
//     eos.RPY2Cs(rhog, p_int, y, cs_ext);
//     amrex::Real eg;
//     eos.RYP2E(rhog, y, p_ext, eg);

//     // amrex::Real p_int;
//     // eos.RTY2P(s_int[URHO], T, y, p_int);
//     // amrex::Real rhog = s_int[URHO] * p_ext / p_int;
//     // amrex::Real eg;
//     // eos.RYP2E(rhog, y, p_ext, eg);

//     // eos.PYT2E(p_ext, y, T, eg); // set p_ghost = p_ext

//     s_ext[URHO] = rhog;
//     s_ext[UMX] = rhog * rhoinv * (dir == 0 ? std::abs(s_int[UMX]) : s_int[UMX]);
//     s_ext[UMY] = rhog * rhoinv * (dir == 1 ? std::abs(s_int[UMY]) : s_int[UMY]);
//     s_ext[UMZ] = rhog * rhoinv * (dir == 2 ? std::abs(s_int[UMZ]) : s_int[UMZ]);
//     s_ext[UEDEN] = rhog * (eg + ke);
//     for (int n = 0; n < NUM_SPECIES; ++n) s_ext[UFS + n] = rhog * y[n];
//   }
// }

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