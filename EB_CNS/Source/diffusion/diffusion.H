#ifndef CNS_DIFFUSION_K_H_
#define CNS_DIFFUSION_K_H_

#include <AMReX_CONSTANTS.H>
#include <AMReX_FArrayBox.H>

#include <cmath>

#include "diff_ops.H"
#include "index_macros.H"

/// @brief Diffusion fluxes in idir-direction for species and enthalpy.
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
cns_diff_species(amrex::IntVect const& iv, const int idir,
                 amrex::Array4<amrex::Real const> const& q,
                 amrex::Array4<amrex::Real const> const& coeffs,
                 amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& dxinv,
                 amrex::Real flx[NVAR])
{
#if (NUM_SPECIES > 1) // no point doing diffusion for single species
  using amrex::Real;

  const amrex::IntVect ivm(iv - amrex::IntVect::TheDimensionVector(idir));

  // Get massfrac, molefrac, enthalpy
  Real mass1[NUM_SPECIES], mass2[NUM_SPECIES], mole1[NUM_SPECIES], mole2[NUM_SPECIES];
  // Real mass3[NUM_SPECIES], mass4[NUM_SPECIES], mole3[NUM_SPECIES], mole4[NUM_SPECIES];
  Real hi1[NUM_SPECIES], hi2[NUM_SPECIES];
  for (int n = 0; n < NUM_SPECIES; ++n) {
    mass1[n] = q(iv, QFS + n);
    mass2[n] = q(ivm, QFS + n);
    // mass1[n] = q(2 * iv - ivm, QFS + n);
    // mass2[n] = q(iv, QFS + n);
    // mass3[n] = q(ivm, QFS + n);
    // mass4[n] = q(2 * ivm - iv, QFS + n);
  }
  auto eos = pele::physics::PhysicsType::eos();
  eos.Y2X(mass1, mole1);
  eos.Y2X(mass2, mole2);
  // eos.Y2X(mass3, mole3);
  // eos.Y2X(mass4, mole4);

  // Compute species and enthalpy fluxes for ideal EOS
  // Get species/enthalpy diffusion, compute correction vel
  Real T = q(iv, QTEMP);
  eos.RTY2Hi(q(iv, QRHO), T, mass1, hi1);
  T = q(ivm, QTEMP);
  eos.RTY2Hi(q(ivm, QRHO), T, mass2, hi2);

  // constexpr int order = 4;
  // const Real w1 = -1.0 / 24.0;
  // const Real w2 = 27.0 / 24.0;
  // const Real w3 = -27.0 / 24.0;
  // const Real w4 = 1.0 / 24.0;

  Real Vc = 0.0;
  const Real dpdx = dxinv[idir] * (q(iv, QPRES) - q(ivm, QPRES));
  // const Real dpdx = normal_diff<order>(iv, idir, QPRES, q, dxinv);
  const Real dlnp = dpdx / (0.5 * (q(iv, QPRES) + q(ivm, QPRES)));
  for (int n = 0; n < NUM_SPECIES; ++n) {
    const Real Xface = 0.5 * (mole1[n] + mole2[n]);
    const Real Yface = 0.5 * (mass1[n] + mass2[n]);
    // const Real Xface = (-mole1[n] + 9 * mole2[n] + 9 * mole3[n] - mole4[n]) / 16.0;
    // const Real Yface = (-mass1[n] + 9 * mass2[n] + 9 * mass3[n] - mass4[n]) / 16.0;
    const Real hface = 0.5 * (hi1[n] + hi2[n]);
    const Real dXdx = (mole1[n] - mole2[n]) * dxinv[idir]; // 2nd order
    // const Real dXdx =
    //   (w1 * mole1[n] + w2 * mole2[n] + w3 * mole3[n] + w4 * mole4[n]) *
    //   dxinv[idir]; // 4th order
    const Real rhoD_n = 0.5 * (coeffs(iv, CRHOD + n) + coeffs(ivm, CRHOD + n));
    const Real Vd = -rhoD_n * (dXdx + (Xface - Yface) * dlnp);
    // const Real Vd =
    //   -(interp<order>(iv, 0, CRHOD + n, coeffs)) * (dXdx + (Xface - Yface) * dlnp);
    Vc += Vd;
    flx[UFS + n] += Vd;
    flx[UEDEN] += Vd * hface;
  }
  // Add correction velocity to fluxes so sum(Vd) = 0
  for (int n = 0; n < NUM_SPECIES; ++n) {
    const Real Yface = 0.5 * (mass1[n] + mass2[n]);
    // const Real Yface = (-mass1[n] + 9 * mass2[n] + 9 * mass3[n] - mass4[n]) / 16.0;
    const Real hface = 0.5 * (hi1[n] + hi2[n]);
    flx[UFS + n] -= Yface * Vc;
    flx[UEDEN] -= Yface * hface * Vc;
  }
#endif
}

/**
 * @brief Compute diffusion fluxes.
 *
 * @param i,j,k  x, y, z index.
 * @param dir    direction, 0:x, 1:y, 2:z.
 * @param q      primitive variables.
 * @param coeffs transport coefficients.
 * @param dxinv  1/dx.
 * @param[out] flx  output diffusion fluxes.
 */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
cns_diff(const amrex::IntVect& iv, const int d1,
         amrex::Array4<const amrex::Real> const& q,
         amrex::Array4<const amrex::Real> const& coeffs,
         amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& dxinv,
         amrex::Real flx[NVAR]) noexcept
{
  using amrex::Real;
  const amrex::IntVect ivm = iv - amrex::IntVect::TheDimensionVector(d1);

  const int d2 = d1 == 0 ? 1 : 0;
  const int d3 = d1 == 2 ? 1 : 2;
  AMREX_D_TERM(const int QU1 = QU + d1;, const int QU2 = QU + d2;
               , const int QU3 = QU + d3;)
  AMREX_D_TERM(const int UM1 = UMX + d1;, const int UM2 = UMX + d2;
               , const int UM3 = UMX + d3;)

  constexpr int order = 2;

  // Aij = dA_i/dx_j
  const Real dTdn = normal_diff<order>(iv, d1, QTEMP, q, dxinv);
  const Real u11 = normal_diff<order>(iv, d1, QU1, q, dxinv);
#if (AMREX_SPACEDIM >= 2)
  const Real u21 = normal_diff<order>(iv, d1, QU2, q, dxinv);
  const Real u12 = tangent_diff<order>(iv, d1, d2, QU1, q, dxinv);
  const Real u22 = tangent_diff<order>(iv, d1, d2, QU2, q, dxinv);
#endif
#if (AMREX_SPACEDIM == 3)
  const Real u31 = normal_diff<order>(iv, d1, QU3, q, dxinv);
  const Real u13 = tangent_diff<order>(iv, d1, d3, QU1, q, dxinv);
  const Real u33 = tangent_diff<order>(iv, d1, d3, QU3, q, dxinv);
#endif
  const Real divu = AMREX_D_TERM(u11, +u22, +u33);
  const Real muf = interp<order>(iv, d1, CMU, coeffs);
  const Real xif = interp<order>(iv, d1, CXI, coeffs);
  AMREX_D_TERM(Real tau11 = muf * (2.0 * u11 - (2.0 / 3.0) * divu) + xif * divu;
               , Real tau12 = muf * (u12 + u21);, Real tau13 = muf * (u13 + u31);)

  AMREX_D_TERM(flx[UM1] -= tau11;, flx[UM2] -= tau12;, flx[UM3] -= tau13;)
  flx[UEDEN] -= 0.5 * (AMREX_D_TERM((q(iv, QU1) + q(ivm, QU1)) * tau11,
                                        +(q(iv, QU2) + q(ivm, QU2)) * tau12,
                                        +(q(iv, QU3) + q(ivm, QU3)) * tau13) +
                           (coeffs(iv, CLAM) + coeffs(ivm, CLAM)) * dTdn);

  // Species transport
  cns_diff_species(iv, d1, q, coeffs, dxinv, flx);
}

#endif