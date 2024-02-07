#ifndef CNS_CENTRAL_SCHEME_H_
#define CNS_CENTRAL_SCHEME_H_

#include <AMReX_Array4.H>

#include "index_macros.H"

// Second-order Kinetic Energy & Entropy Preserving (KEEP) scheme
// Ref: Fujiwara et al. (2023), https://doi.org/10.1016/j.jcp.2023.111973
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
cns_KEEP2(int i, int j, int k, int dir, amrex::Array4<amrex::Real> const& flx,
          amrex::Array4<const amrex::Real> const& q) noexcept
{
  using amrex::Real;
  const amrex::IntVect iv(AMREX_D_DECL(i, j, k));
  const amrex::IntVect ivm = iv - amrex::IntVect::TheDimensionVector(dir);

  const Real phi = q(iv, QWBAR) / q(ivm, QWBAR) * q(ivm, QRHO) / q(iv, QRHO);

  // Mass: rho * Y * u_i
  flx(iv, URHO) = 0.0;
  for (int n = 0; n < NUM_SPECIES; ++n) {
    flx(iv, UFS + n) =
      0.25 *
      (q(ivm, QRHO) * q(ivm, QFS + n) / phi + q(iv, QRHO) * q(iv, QFS + n) * phi) *
      (q(iv, QU + dir) + q(ivm, QU + dir));
    flx(iv, URHO) += flx(iv, UFS + n);
  }

  // Momentum: rho * u_i * u_j + p * delta_ij
  flx(iv, UMX) = 0.125 * (q(ivm, QRHO) / phi + q(iv, QRHO) * phi) *
                 (q(ivm, QU + dir) + q(iv, QU + dir)) * (q(ivm, QU) + q(iv, QU));
#if AMREX_SPACEDIM >= 2
  flx(iv, UMY) = 0.125 * (q(ivm, QRHO) / phi + q(iv, QRHO) * phi) *
                 (q(ivm, QU + dir) + q(iv, QU + dir)) * (q(ivm, QV) + q(iv, QV));
#endif
#if AMREX_SPACEDIM == 3
  flx(iv, UMZ) = 0.125 * (q(ivm, QRHO) / phi + q(iv, QRHO) * phi) *
                 (q(ivm, QU + dir) + q(iv, QU + dir)) * (q(ivm, QW) + q(iv, QW));
#endif
  flx(iv, UMX + dir) += 0.5 * (q(ivm, QPRES) + q(iv, QPRES));

  // Energy: rho * (ke + ie) * u_i + p * u_i
  flx(iv, UEDEN) = 0.125 * (q(ivm, QRHO) / phi + q(iv, QRHO) * phi) *
                     (q(ivm, QU + dir) + q(iv, QU + dir)) *
                     (AMREX_D_TERM(q(ivm, QU) * q(iv, QU), +q(ivm, QV) * q(iv, QV),
                                   +q(ivm, QW) * q(iv, QW))) // ke
                   + 0.25 *
                       (q(ivm, QRHO) * q(ivm, QEINT) + q(iv, QRHO) * q(iv, QEINT)) *
                       (q(ivm, QU + dir) + q(iv, QU + dir)) // ie
                   + 0.5 * (q(iv, QPRES) * q(ivm, QU + dir) +
                            q(ivm, QPRES) * q(iv, QU + dir)); // p * u_i
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
cns_KEEP_base(amrex::Real rhom, amrex::Real rho,
              amrex::Real Ym[NUM_SPECIES], amrex::Real Y[NUM_SPECIES],
              amrex::Real um, amrex::Real u, amrex::Real vm, amrex::Real v, 
              amrex::Real wm, amrex::Real w, amrex::Real pm, amrex::Real p,
              amrex::Real eim, amrex::Real ei, amrex::Real wbarm, amrex::Real wbar, 
              amrex::Real flx[NVAR]) noexcept
{
  const amrex::Real phi = wbar / wbarm * rhom / rho;

  flx[URHO] = 0.25 * (rhom / phi + rho * phi) * (um + u);
  for (int n = 0; n < NUM_SPECIES; ++n) {
    flx[UFS + n] = 0.25 * (rhom * Ym[n] / phi + rho * Y[n] * phi) * (um + u);
  }
  flx[UMX] = 0.5 * flx[URHO] * (um + u) + 0.5 * (pm + p);
#if AMREX_SPACEDIM >= 2
  flx[UMY] = 0.5 * flx[URHO] * (vm + v);
#endif
#if AMREX_SPACEDIM == 3
  flx[UMZ] = 0.5 * flx[URHO] * (wm + w);
#endif
  flx[UEDEN] = 0.5 * flx[URHO] * (AMREX_D_TERM(um * u, +vm * v, +wm * w)) // ke
               + 0.25 * (rhom * eim + rho * ei) * (um + u)                // ie
               + 0.5 * (p * um + pm * u);                                 // p * u_i
}

// Fourth-order Kinetic Energy & Entropy Preserving (KEEP) scheme
// Ref: Kuya & Kawai (2021), http://doi.org/10.1016/j.jcp.2021.110482
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
cns_KEEP4(int i, int j, int k, int dir, amrex::Array4<amrex::Real> const& flx,
          amrex::Array4<const amrex::Real> const& q) noexcept
{
  using amrex::Real;
  const amrex::IntVect iv(AMREX_D_DECL(i, j, k));

  const int QUN = QU + dir;
  const int QUT = QU + (dir + 1) % AMREX_SPACEDIM;
  const int QUTT = QU + (dir + 2) % AMREX_SPACEDIM;
  
  // (KEEP4 flux) = 4/3 * (KEEP2 flux)_{i-1,i} - 1/6 * ((KEEP2 flux)_{i-1,i+1} + (KEEP2 flux)_{i-2,i})
  Real flx_tmp[NVAR], Ym[NUM_SPECIES], Y[NUM_SPECIES];
  
  // (KEEP2 flux)_{i-1,i}
  const amrex::IntVect ivm = iv - amrex::IntVect::TheDimensionVector(dir);
  for (int n = 0; n < NUM_SPECIES; ++n) {
    Ym[n] = q(ivm, QFS + n);
    Y[n] = q(iv, QFS + n);
  }
  cns_KEEP_base(q(ivm, QRHO), q(iv, QRHO), Ym, Y, q(ivm, QUN), q(iv, QUN),
                q(ivm, QUT), q(iv, QUT), q(ivm, QUTT), q(iv, QUTT), q(ivm, QPRES),
                q(iv, QPRES), q(ivm, QEINT), q(iv, QEINT), q(ivm, QWBAR),
                q(iv, QWBAR), flx_tmp);
  for (int n = 0; n < NUM_SPECIES; ++n) {
    flx(iv, UFS + n) = flx_tmp[UFS + n] * 4.0 / 3.0;
  }
  flx(iv, URHO) = flx_tmp[URHO] * 4.0 / 3.0;
  AMREX_D_TERM(flx(iv, QUN) = flx_tmp[UMX] * 4.0 / 3.0; // we can do this because UMX = QU
               , flx(iv, QUT) = flx_tmp[UMY] * 4.0 / 3.0;
               , flx(iv, QUTT) = flx_tmp[UMZ] * 4.0 / 3.0;)
  flx(iv, UEDEN) = flx_tmp[UEDEN] * 4.0 / 3.0;

  // (KEEP2 flux)_{i-1,i+1}
  const amrex::IntVect ivp = iv + amrex::IntVect::TheDimensionVector(dir);
  for (int n = 0; n < NUM_SPECIES; ++n) {
    // Ym[n] = q(ivm, QFS + n);
    Y[n] = q(ivp, QFS + n);
  }
  cns_KEEP_base(q(ivm, QRHO), q(ivp, QRHO), Ym, Y, q(ivm, QUN), q(ivp, QUN),
                q(ivm, QUT), q(ivp, QUT), q(ivm, QUTT), q(ivp, QUTT), q(ivm, QPRES),
                q(ivp, QPRES), q(ivm, QEINT), q(ivp, QEINT), q(ivm, QWBAR),
                q(ivp, QWBAR), flx_tmp);
  for (int n = 0; n < NUM_SPECIES; ++n) {
    flx(iv, UFS + n) -= flx_tmp[UFS + n] / 6.0;
  }
  flx(iv, URHO) -= flx_tmp[URHO] / 6.0;
  AMREX_D_TERM(flx(iv, QUN) -= flx_tmp[UMX] / 6.0;
               , flx(iv, QUT) -= flx_tmp[UMY] / 6.0;
               , flx(iv, QUTT) -= flx_tmp[UMZ] / 6.0;)
  flx(iv, UEDEN) -= flx_tmp[UEDEN] / 6.0;

  // (KEEP2 flux)_{i-2,i}
  const amrex::IntVect ivmm = iv - 2 * amrex::IntVect::TheDimensionVector(dir);
  for (int n = 0; n < NUM_SPECIES; ++n) {
    Ym[n] = q(ivmm, QFS + n);
    Y[n] = q(iv, QFS + n);
  }
  cns_KEEP_base(q(ivmm, QRHO), q(iv, QRHO), Ym, Y, q(ivmm, QUN), q(iv, QUN),
                q(ivmm, QUT), q(iv, QUT), q(ivmm, QUTT), q(iv, QUTT), q(ivmm, QPRES),
                q(iv, QPRES), q(ivmm, QEINT), q(iv, QEINT), q(ivmm, QWBAR),
                q(iv, QWBAR), flx_tmp);
  for (int n = 0; n < NUM_SPECIES; ++n) {
    flx(iv, UFS + n) -= flx_tmp[UFS + n] / 6.0;
  }
  flx(iv, URHO) -= flx_tmp[URHO] / 6.0;
  AMREX_D_TERM(flx(iv, QUN) -= flx_tmp[UMX] / 6.0;
               , flx(iv, QUT) -= flx_tmp[UMY] / 6.0;
               , flx(iv, QUTT) -= flx_tmp[UMZ] / 6.0;)
  flx(iv, UEDEN) -= flx_tmp[UEDEN] / 6.0;
}

// // Fourth-order Kinetic Energy & Entropy Preserving (KEEP) scheme
// AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
// cns_KEEP4(int i, int j, int k, int dir, amrex::Array4<amrex::Real> const& flx,
//           amrex::Array4<const amrex::Real> const& q) noexcept
// {
//   using amrex::Real;
//   amrex::IntVect iv(AMREX_D_DECL(i, j, k));
//   amrex::IntVect ivm = iv - amrex::IntVect::TheDimensionVector(dir);

//   // Get KEEP2 fluxes at i-1/2, store in flx
//   Real phi = q(iv, QWBAR) / q(ivm, QWBAR) * q(ivm, QRHO) / q(iv, QRHO);

//   // Mass: rho * Y * u_i
//   flx(iv, URHO) = 0.0;
//   for (int n = 0; n < NUM_SPECIES; ++n) {
//     flx(iv, UFS + n) =
//       0.25 *
//       (q(ivm, QRHO) * q(ivm, QFS + n) / phi + q(iv, QRHO) * q(iv, QFS + n) * phi) *
//       (q(iv, QU + dir) + q(ivm, QU + dir));
//     flx(iv, URHO) += flx(iv, UFS + n);
//   }

//   // Momentum: rho * u_i * u_j + p * delta_ij
//   flx(iv, UMX) = 0.125 * (q(ivm, QRHO) / phi + q(iv, QRHO) * phi) *
//                  (q(ivm, QU + dir) + q(iv, QU + dir)) * (q(ivm, QU) + q(iv, QU));
// #if AMREX_SPACEDIM >= 2
//   flx(iv, UMY) = 0.125 * (q(ivm, QRHO) / phi + q(iv, QRHO) * phi) *
//                  (q(ivm, QU + dir) + q(iv, QU + dir)) * (q(ivm, QV) + q(iv, QV));
// #endif
// #if AMREX_SPACEDIM == 3
//   flx(iv, UMZ) = 0.125 * (q(ivm, QRHO) / phi + q(iv, QRHO) * phi) *
//                  (q(ivm, QU + dir) + q(iv, QU + dir)) * (q(ivm, QW) + q(iv, QW));
// #endif
//   flx(iv, UMX + dir) += 0.5 * (q(ivm, QPRES) + q(iv, QPRES));

//   // Energy: rho * (ke + ie) * u_i + p * u_i
//   flx(iv, UEDEN) = 0.125 * (q(ivm, QRHO) / phi + q(iv, QRHO) * phi) *
//                      (q(ivm, QU + dir) + q(iv, QU + dir)) *
//                      (AMREX_D_TERM(q(ivm, QU) * q(iv, QU), +q(ivm, QV) * q(iv, QV),
//                                    +q(ivm, QW) * q(iv, QW))) // ke
//                    + 0.25 *
//                        (q(ivm, QRHO) * q(ivm, QEINT) + q(iv, QRHO) * q(iv, QEINT)) *
//                        (q(ivm, QU + dir) + q(iv, QU + dir)) // ie
//                    + 0.5 * (q(iv, QPRES) * q(ivm, QU + dir) +
//                             q(ivm, QPRES) * q(iv, QU + dir)); // p * u_i

//   // Get KEEP2 fluxes at i+1/2, store in tmp
//   Real tmp[NVAR];  

//   // shift one cell in dir direction
//   iv += amrex::IntVect::TheDimensionVector(dir);
//   ivm += amrex::IntVect::TheDimensionVector(dir);
  
//   // these are the same as above
//   phi = q(iv, QWBAR) / q(ivm, QWBAR) * q(ivm, QRHO) / q(iv, QRHO);

//   tmp[URHO] = 0.0;
//   for (int n = 0; n < NUM_SPECIES; ++n) {
//     tmp[UFS + n] =
//       0.125 *
//       (q(ivm, QRHO) * q(ivm, QFS + n) / phi + q(iv, QRHO) * q(iv, QFS + n) * phi) *
//       (q(iv, QU + dir) + q(ivm, QU + dir));
//     tmp[URHO] += tmp[UFS + n];
//   }

//   tmp[UMX] = 0.125 * (q(ivm, QRHO) / phi + q(iv, QRHO) * phi) *
//              (q(ivm, QU + dir) + q(iv, QU + dir)) * (q(ivm, QU) + q(iv, QU));
// #if AMREX_SPACEDIM >= 2
//   tmp[UMY] = 0.125 * (q(ivm, QRHO) / phi + q(iv, QRHO) * phi) *
//              (q(ivm, QU + dir) + q(iv, QU + dir)) * (q(ivm, QV) + q(iv, QV));
// #endif
// #if AMREX_SPACEDIM == 3
//   tmp[UMZ] = 0.125 * (q(ivm, QRHO) / phi + q(iv, QRHO) * phi) *
//              (q(ivm, QU + dir) + q(iv, QU + dir)) * (q(ivm, QW) + q(iv, QW));
// #endif

//   tmp[UEDEN] =
//     0.0125 * (q(ivm, QRHO) / phi + q(iv, QRHO) * phi) *
//       (q(ivm, QU + dir) + q(iv, QU + dir)) *
//       (AMREX_D_TERM(q(ivm, QU) * q(iv, QU), +q(ivm, QV) * q(iv, QV),
//                     +q(ivm, QW) * q(iv, QW))) +
//     0.25 * (q(ivm, QRHO) * q(ivm, QEINT) + q(iv, QRHO) * q(iv, QEINT)) *
//       (q(ivm, QU + dir) + q(iv, QU + dir)) +
//     0.5 * (q(iv, QPRES) * q(ivm, QU + dir) + q(ivm, QPRES) * q(iv, QU + dir));

//   // Get KEEP4 fluxes: flx_{i-0.5} = 4/3 * flx_{i-0.5} - 1/6 * flx_{i+0.5}
//   iv = amrex::IntVect(AMREX_D_DECL(i, j, k));  
//   for (int n = 0; n < NVAR; ++n) {
//     flx(iv, n) = flx(iv, n) * 4.0 / 3.0 - tmp[n] / 6.0;
//   }
// }

// Conservative filtering
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
cns_consFilter(int i, int j, int k, amrex::Array4<amrex::Real> const& flx,
               amrex::Array4<amrex::Real> const& dSdt,
               amrex::Array4<const amrex::Real> const& s) noexcept
{
  using amrex::Real;
  const amrex::IntVect iv(AMREX_D_DECL(i, j, k));
  const amrex::IntVect ivm = iv - amrex::IntVect::TheDimensionVector(0);

  // Compute S^{n+1}

  // Modify fluxes

  // Update dSdt (OR need to wait until all dir are done?)

}

#endif