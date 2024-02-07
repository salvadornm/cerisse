#ifndef CNS_DIFFUSION_K_H_
#define CNS_DIFFUSION_K_H_

#include <AMReX_CONSTANTS.H>
#include <AMReX_FArrayBox.H>

#include <cmath>

#include "index_macros.H"
#include "diff_ops.H"

/// @brief Diffusion fluxes in idir-direction for species and enthalpy.
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
cns_diff_species(amrex::IntVect const& iv, const int idir, amrex::Array4<amrex::Real const> const& q,
                 amrex::Array4<amrex::Real const> const& coeffs,
                 amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& dxinv,
                 amrex::Array4<amrex::Real> const& flx) 
{
#if (NUM_SPECIES > 1) // no point doing diffusion for single species
  using amrex::Real;
  
  const amrex::IntVect ivm(iv - amrex::IntVect::TheDimensionVector(idir));

  // Get massfrac, molefrac, enthalpy
  Real mass1[NUM_SPECIES], mass2[NUM_SPECIES], mass3[NUM_SPECIES], mass4[NUM_SPECIES];
  Real mole1[NUM_SPECIES], mole2[NUM_SPECIES], mole3[NUM_SPECIES], mole4[NUM_SPECIES];
  Real hi1[NUM_SPECIES], hi2[NUM_SPECIES];
  for (int n = 0; n < NUM_SPECIES; ++n) {
    // mass1[n] = q(iv, QFS + n);
    // mass2[n] = q(ivm, QFS + n);
    mass1[n] = q(2 * iv - ivm, QFS + n);
    mass2[n] = q(iv, QFS + n);
    mass3[n] = q(ivm, QFS + n);
    mass4[n] = q(2 * ivm - iv, QFS + n);
  }
  auto eos = pele::physics::PhysicsType::eos();
  eos.Y2X(mass1, mole1);
  eos.Y2X(mass2, mole2);
  eos.Y2X(mass3, mole3);
  eos.Y2X(mass4, mole4);

  // Compute species and enthalpy fluxes for ideal EOS
  // Get species/enthalpy diffusion, compute correction vel
  Real T = q(iv, QTEMP);
  eos.T2Hi(T, hi1);
  T = q(ivm, QTEMP);
  eos.T2Hi(T, hi2);

  constexpr int order = 4;
  const Real w1 = -1.0 / 24.0;
  const Real w2 = 27.0 / 24.0;
  const Real w3 = -27.0 / 24.0;
  const Real w4 = 1.0 / 24.0;

  Real Vc = 0.0;
  const Real dpdx = normal_diff<order>(iv, idir, QPRES, q, dxinv);
  const Real dlnp = dpdx / (0.5 * (q(iv, QPRES) + q(ivm, QPRES)));
  for (int n = 0; n < NUM_SPECIES; ++n) {
    // const Real Xface = 0.5 * (mole1[n] + mole2[n]);
    // const Real Yface = 0.5 * (mass1[n] + mass2[n]);
    const Real Xface = (-mole1[n] + 9 * mole2[n] + 9 * mole3[n] - mole4[n]) / 16.0;
    const Real Yface = (-mass1[n] + 9 * mass2[n] + 9 * mass3[n] - mass4[n]) / 16.0;
    const Real hface = 0.5 * (hi1[n] + hi2[n]);
    // const Real dXdx = (mole1[n] - mole2[n]) * dxinv[idir]; // 2nd order
    const Real dXdx =
      (w1 * mole1[n] + w2 * mole2[n] + w3 * mole3[n] + w4 * mole4[n]) *
      dxinv[idir]; // 4th order
    const Real Vd = -(interp<order>(iv, 0, CRHOD + n, coeffs)) *
                           (dXdx + (Xface - Yface) * dlnp);
    Vc += Vd;
    flx(iv, UFS + n) += Vd;
    flx(iv, UEDEN) += Vd * hface;
  }
  // Add correction velocity to fluxes so sum(Vd) = 0
  for (int n = 0; n < NUM_SPECIES; ++n) {
    // const Real Yface = 0.5 * (mass1[n] + mass2[n]);
    const Real Yface = (-mass1[n] + 9 * mass2[n] + 9 * mass3[n] - mass4[n]) / 16.0;
    const Real hface = 0.5 * (hi1[n] + hi2[n]);
    flx(iv, UFS + n) -= Yface * Vc;
    flx(iv, UEDEN) -= Yface * hface * Vc;
  }
#endif
}

/// @brief Diffusion fluxes in x-direction. 4th order central difference.
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
cns_diff_x(int i, int j, int k, amrex::Array4<amrex::Real const> const& q,
           amrex::Array4<amrex::Real const> const& coeffs,
           amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& dxinv,
           amrex::Array4<amrex::Real> const& fx) noexcept
{
  using amrex::Real;
  const amrex::IntVect iv{AMREX_D_DECL(i, j, k)};
  constexpr int order = 4;

  Real dTdx = normal_diff<order>(iv, 0, QTEMP, q, dxinv);
  AMREX_D_TERM(Real dudx = normal_diff<order>(iv, 0, QU, q, dxinv);
               , Real dvdx = normal_diff<order>(iv, 0, QV, q, dxinv);
               , Real dwdx = normal_diff<order>(iv, 0, QW, q, dxinv);)
#if (AMREX_SPACEDIM >= 2)
  Real dudy = tangent_diff<order>(iv, 0, 1, QU, q, dxinv);
  Real dvdy = tangent_diff<order>(iv, 0, 1, QV, q, dxinv);
#endif
#if (AMREX_SPACEDIM == 3)
  Real dudz = tangent_diff<order>(iv, 0, 2, QU, q, dxinv);
  Real dwdz = tangent_diff<order>(iv, 0, 2, QW, q, dxinv);
#endif
  Real divu = AMREX_D_TERM(dudx, +dvdy, +dwdz);

  const Real mu = interp<order>(iv, 0, CMU, coeffs);
  const Real xi = interp<order>(iv, 0, CXI, coeffs);

  AMREX_D_TERM(Real tauxx = mu * (2.0 * dudx - 2.0 / 3.0 * divu) + xi * divu;
               , Real tauxy = mu * (dudy + dvdx);, Real tauxz = mu * (dudz + dwdx););

  AMREX_D_TERM(fx(i, j, k, UMX) += -tauxx;, fx(i, j, k, UMY) += -tauxy;
               , fx(i, j, k, UMZ) += -tauxz;);
  fx(i, j, k, UEDEN) += -(AMREX_D_TERM(interp<order>(iv, 0, QU, q) * tauxx,
                                       +interp<order>(iv, 0, QV, q) * tauxy,
                                       +interp<order>(iv, 0, QW, q) * tauxz)) -
                        interp<order>(iv, 0, CLAM, coeffs) * dTdx;

  // Species transport
  cns_diff_species(iv, 0, q, coeffs, dxinv, fx);
}

#if (AMREX_SPACEDIM >= 2)
/// @brief Diffusion fluxes in y-direction. 4th order central difference.
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
cns_diff_y(int i, int j, int k, amrex::Array4<amrex::Real const> const& q,
           amrex::Array4<amrex::Real const> const& coeffs,
           amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& dxinv,
           amrex::Array4<amrex::Real> const& fy) noexcept
{
  using amrex::Real;
  const amrex::IntVect iv{AMREX_D_DECL(i, j, k)};
  constexpr int order = 4;

  Real dTdy = normal_diff<order>(iv, 1, QTEMP, q, dxinv);
  AMREX_D_TERM(Real dudy = normal_diff<order>(iv, 1, QU, q, dxinv);
               , Real dvdy = normal_diff<order>(iv, 1, QV, q, dxinv);
               , Real dwdy = normal_diff<order>(iv, 1, QW, q, dxinv);)

  Real dudx = tangent_diff<order>(iv, 1, 0, QU, q, dxinv);
  Real dvdx = tangent_diff<order>(iv, 1, 0, QV, q, dxinv);
#if (AMREX_SPACEDIM == 3)
  Real dvdz = tangent_diff<order>(iv, 1, 2, QV, q, dxinv);
  Real dwdz = tangent_diff<order>(iv, 1, 2, QW, q, dxinv);
#endif
  Real divu = AMREX_D_TERM(dudx, +dvdy, +dwdz);

  const Real mu = interp<order>(iv, 1, CMU, coeffs);
  const Real xi = interp<order>(iv, 1, CXI, coeffs);

  AMREX_D_TERM(Real tauyy = mu * (2.0 * dvdy - (2.0 / 3.0) * divu) + xi * divu;
               , Real tauxy = mu * (dudy + dvdx);, Real tauyz = mu * (dwdy + dvdz););

  AMREX_D_TERM(fy(i, j, k, UMX) += -tauxy;, fy(i, j, k, UMY) += -tauyy;
               , fy(i, j, k, UMZ) += -tauyz;);
  fy(i, j, k, UEDEN) += -(AMREX_D_TERM(interp<order>(iv, 1, QU, q) * tauxy,
                                       +interp<order>(iv, 1, QV, q) * tauyy,
                                       +interp<order>(iv, 1, QW, q) * tauyz)) -
                        interp<order>(iv, 1, CLAM, coeffs) * dTdy;

  // Species transport
  cns_diff_species(iv, 1, q, coeffs, dxinv, fy);
}
#endif

#if (AMREX_SPACEDIM == 3)
/// @brief Diffusion fluxes in z-direction. 4th order central difference.
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
cns_diff_z(int i, int j, int k, amrex::Array4<amrex::Real const> const& q,
           amrex::Array4<amrex::Real const> const& coeffs,
           amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& dxinv,
           amrex::Array4<amrex::Real> const& fz) noexcept
{
  using amrex::Real;
  const amrex::IntVect iv{AMREX_D_DECL(i, j, k)};
  constexpr int order = 4;

  Real dTdz = normal_diff<order>(iv, 2, QTEMP, q, dxinv);
  Real dudz = normal_diff<order>(iv, 2, QU, q, dxinv);
  Real dvdz = normal_diff<order>(iv, 2, QV, q, dxinv);
  Real dwdz = normal_diff<order>(iv, 2, QW, q, dxinv);

  Real dudx = tangent_diff<order>(iv, 2, 0, QU, q, dxinv);
  Real dwdx = tangent_diff<order>(iv, 2, 0, QW, q, dxinv);
  Real dvdy = tangent_diff<order>(iv, 2, 1, QV, q, dxinv);
  Real dwdy = tangent_diff<order>(iv, 2, 1, QW, q, dxinv);

  Real divu = dudx + dvdy + dwdz;
  const Real mu = interp<order>(iv, 2, CMU, coeffs);
  const Real xi = interp<order>(iv, 2, CXI, coeffs);
  
  Real tauxz = mu * (dudz + dwdx);
  Real tauyz = mu * (dvdz + dwdy);
  Real tauzz = mu * (2.0 * dwdz - (2.0 / 3.0) * divu) + xi * divu;

  fz(i, j, k, UMX) += -tauxz;
  fz(i, j, k, UMY) += -tauyz;
  fz(i, j, k, UMZ) += -tauzz;
  fz(i, j, k, UEDEN) +=
    -(interp<order>(iv, 2, QU, q) * tauxz + interp<order>(iv, 2, QV, q) * tauyz +
      interp<order>(iv, 2, QW, q) * tauzz) -
    interp<order>(iv, 2, CLAM, coeffs) * dTdz;

  // Species transport
  cns_diff_species(iv, 2, q, coeffs, dxinv, fz);
}
#endif

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
cns_diff(int i, int j, int k, int dir, amrex::Array4<const amrex::Real> const& q,
         amrex::Array4<const amrex::Real> const& coeffs,
         amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& dxinv,
         amrex::Array4<amrex::Real> const& flx) noexcept
{
  if (dir == 0) {
    cns_diff_x(i, j, k, q, coeffs, dxinv, flx);
  }
#if (AMREX_SPACEDIM >= 2)
  else if (dir == 1) {
    cns_diff_y(i, j, k, q, coeffs, dxinv, flx);
  }
#endif
#if (AMREX_SPACEDIM == 3)
  else {
    cns_diff_z(i, j, k, q, coeffs, dxinv, flx);
  }
#endif
}

#endif