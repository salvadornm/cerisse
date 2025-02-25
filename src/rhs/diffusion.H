#ifndef CNS_DIFFUSION_K_H_
#define CNS_DIFFUSION_K_H_

#include <AMReX_CONSTANTS.H>
#include <AMReX_FArrayBox.H>

#include <cmath>
#include "diff_ops.H"


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

#endif