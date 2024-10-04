#ifndef DiffusionCD_H_
#define DiffusionCD_H_

#include <AMReX_FArrayBox.H>

/////////////////////////////////////////////////////////////////////////////
// Normal derivative with variable order using central difference.
// order: order of accuracy.
// iv: cell index (i, j, k).
// idir: direction of derivative.
// comp: component of the field.
// q: primitive variable Array4.
// dxinv: inverse of cell size.
template <int order>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real normal_diff(
    amrex::IntVect iv, int idir, int comp,
    amrex::Array4<amrex::Real const> const& q,
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& dxinv) noexcept;

/// @brief 2nd order central difference.
template <>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real normal_diff<2>(
    amrex::IntVect iv, int idir, int comp,
    amrex::Array4<amrex::Real const> const& q,
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& dxinv) noexcept {
  const auto ivn = amrex::IntVect::TheDimensionVector(idir);
  return (q(iv, comp) - q(iv - ivn, comp)) * dxinv[idir];
}

// Tangential derivative with variable order using central difference.
template <int order>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real tangent_diff(
    amrex::IntVect iv, int idir, int tdir, int comp,
    amrex::Array4<amrex::Real const> const& q,
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& dxinv) noexcept;

/// @brief 2nd order central difference. d/d(tdir) on face normal to idir.
template <>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real tangent_diff<2>(
    amrex::IntVect iv, int idir, int tdir, int comp,
    amrex::Array4<amrex::Real const> const& q,
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& dxinv) noexcept {
  const auto ivn = amrex::IntVect::TheDimensionVector(idir);
  const auto ivt = amrex::IntVect::TheDimensionVector(tdir);
  return (q(iv + ivt, comp) - q(iv - ivt, comp) + q(iv - ivn + ivt, comp) -
          q(iv - ivn - ivt, comp)) *
         0.25 * dxinv[tdir];
}

// Interpolation with variable order using polynomial.
template <int order>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real interp(
    amrex::IntVect iv, int idir, int comp,
    amrex::Array4<amrex::Real const> const& q) noexcept;

template <>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real interp<2>(
    amrex::IntVect iv, int idir, int comp,
    amrex::Array4<amrex::Real const> const& q) noexcept {
  const auto ivn = amrex::IntVect::TheDimensionVector(idir);
  return 0.5 * (q(iv, comp) + q(iv - ivn, comp));
}
/////////////////////////////////////////////////////////////////////////////

template <typename cls_t, typename... visc_types>
class diffusion_t : public visc_types... {
 public:
  // AMREX_GPU_HOST_DEVICE
  diffusion_t() {}

  // AMREX_GPU_HOST_DEVICE
  ~diffusion_t() {}

  void inline dflux(const Geometry& geom, const MFIter& mfi,
                    const Array4<Real>& prims, const Array4<Real>& flx,
                    const Array4<Real>& rhs, const cls_t* cls) {
    const GpuArray<Real, AMREX_SPACEDIM> dxinv = geom.InvCellSizeArray();
    const Box& bx = mfi.tilebox();

    // Get coefficients
    const Box& bxg1 = amrex::grow(bx, 1);
    FArrayBox coeffs(bxg1, 3 + cls_t::NSPEC, The_Async_Arena());
    const auto& mu_arr = coeffs.array(0);
    const auto& xi_arr = coeffs.array(1);
    const auto& lam_arr = coeffs.array(2);
    const auto& rD_arr = coeffs.array(3);
    amrex::ParallelFor(
        bxg1, [=, *this] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
          this->visc(i, j, k, mu_arr, prims, *cls);
          this->bulk_visc(i, j, k, xi_arr, prims, *cls);
          this->cond(i, j, k, lam_arr, prims, *cls);
          this->spec_diff(i, j, k, rD_arr, prims, *cls);
        });

    // x-direction
    int cdir = 0;
    const Box& xflxbx = amrex::surroundingNodes(bx, cdir);
    amrex::ParallelFor(
        xflxbx, [=, *this] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
          for (int n = 0; n < cls_t::NCONS; n++) flx(i, j, k, n) = 0.0;
          this->cns_diff_x(i, j, k, flx, prims, mu_arr, xi_arr, lam_arr, rD_arr,
                           dxinv, *cls);
        });
    amrex::ParallelFor(
        bx, cls_t::NCONS,
        [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
          rhs(i, j, k, n) +=
              dxinv[cdir] * (flx(i, j, k, n) - flx(i + 1, j, k, n));
        });

#if AMREX_SPACEDIM >= 2
    // y-direction
    cdir = 1;
    const Box& yflxbx = amrex::surroundingNodes(bx, cdir);
    amrex::ParallelFor(
        yflxbx, [=, *this] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
          for (int n = 0; n < cls_t::NCONS; n++) flx(i, j, k, n) = 0.0;
          this->cns_diff_y(i, j, k, flx, prims, mu_arr, xi_arr, lam_arr, rD_arr,
                           dxinv, *cls);
        });
    amrex::ParallelFor(
        bx, cls_t::NCONS,
        [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
          rhs(i, j, k, n) +=
              dxinv[cdir] * (flx(i, j, k, n) - flx(i, j + 1, k, n));
        });
#endif

#if AMREX_SPACEDIM == 3
    // z-direction
    cdir = 2;
    const Box& zflxbx = amrex::surroundingNodes(bx, cdir);
    amrex::ParallelFor(
        zflxbx, [=, *this] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
          for (int n = 0; n < cls_t::NCONS; n++) flx(i, j, k, n) = 0.0;
          this->cns_diff_z(i, j, k, flx, prims, mu_arr, xi_arr, lam_arr, rD_arr,
                           dxinv, *cls);
        });
    amrex::ParallelFor(
        bx, cls_t::NCONS,
        [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
          rhs(i, j, k, n) +=
              dxinv[cdir] * (flx(i, j, k, n) - flx(i, j, k + 1, n));
        });
#endif
  }

  AMREX_GPU_DEVICE AMREX_FORCE_INLINE void cns_diff_species(
      amrex::IntVect const& iv, const int idir,
      amrex::Array4<amrex::Real const> const& q,
      amrex::Array4<amrex::Real const> const& rD_arr,
      amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& dxinv,
      amrex::Array4<amrex::Real> const& flx, cls_t const& cls) const {
#if (NUM_SPECIES > 1)  // no point doing diffusion for single species
    using amrex::Real;
    const amrex::IntVect ivm(iv - amrex::IntVect::TheDimensionVector(idir));

    // Get massfrac, molefrac, enthalpy
    Real mass1[NUM_SPECIES], mass2[NUM_SPECIES];
    Real mole1[NUM_SPECIES], mole2[NUM_SPECIES];
    Real hi1[NUM_SPECIES], hi2[NUM_SPECIES];
    for (int n = 0; n < NUM_SPECIES; ++n) {
      mass1[n] = q(iv, cls.QFS + n);
      mass2[n] = q(ivm, cls.QFS + n);
    }
    auto eos = pele::physics::PhysicsType::eos();
    eos.Y2X(mass1, mole1);
    eos.Y2X(mass2, mole2);

    // Compute species and enthalpy fluxes for ideal EOS
    // Get species/enthalpy diffusion, compute correction vel
    Real T = q(iv, cls.QT);
    eos.T2Hi(T, hi1);
    T = q(ivm, cls.QT);
    eos.T2Hi(T, hi2);
    // Change to SI units
    for (int n = 0; n < NUM_SPECIES; ++n) {
      hi1[n] *= 1.e-7;
      hi2[n] *= 1.e-7;  // TODO: correct?
    }

    // 2nd order

    Real Vc = 0.0;
    const Real dpdx = (q(iv, cls.QPRES) - q(ivm, cls.QPRES)) * dxinv[idir];
    const Real dlnp = dpdx / (0.5 * (q(iv, cls.QPRES) + q(ivm, cls.QPRES)));
    for (int n = 0; n < NUM_SPECIES; ++n) {
      const Real Xface = 0.5 * (mole1[n] + mole2[n]);
      const Real Yface = 0.5 * (mass1[n] + mass2[n]);
      const Real hface = 0.5 * (hi1[n] + hi2[n]);
      const Real dXdx = (mole1[n] - mole2[n]) * dxinv[idir];
      const Real Vd = -0.5 * (rD_arr(iv, n) + rD_arr(ivm, n)) *
                      (dXdx + (Xface - Yface) * dlnp);
      Vc += Vd;
      flx(iv, cls.UFS + n) += Vd;
      flx(iv, cls.UET) += Vd * hface;
    }
    // Add correction velocity to fluxes so sum(Vd) = 0
    for (int n = 0; n < NUM_SPECIES; ++n) {
      const Real Yface = 0.5 * (mass1[n] + mass2[n]);
      const Real hface = 0.5 * (hi1[n] + hi2[n]);
      flx(iv, cls.UFS + n) -= Yface * Vc;
      flx(iv, cls.UET) -= Yface * hface * Vc;
    }
#endif
  }

  AMREX_GPU_DEVICE AMREX_FORCE_INLINE void cns_diff_x(
      int i, int j, int k, amrex::Array4<amrex::Real> const& fx,
      amrex::Array4<const amrex::Real> const& prims,
      amrex::Array4<const amrex::Real> const& mu_arr,
      amrex::Array4<const amrex::Real> const& xi_arr,
      amrex::Array4<const amrex::Real> const& lam_arr,
      amrex::Array4<const amrex::Real> const& rD_arr,
      amrex::GpuArray<amrex::Real, amrex::SpaceDim> const& dxinv,
      cls_t const& cls) const {
    using amrex::Real;
    const amrex::IntVect iv{AMREX_D_DECL(i, j, k)};
    constexpr int order = 2;

    Real dTdx = normal_diff<order>(iv, 0, cls.QT, prims, dxinv);
    AMREX_D_TERM(Real dudx = normal_diff<order>(iv, 0, cls.QU, prims, dxinv);
                 , Real dvdx = normal_diff<order>(iv, 0, cls.QV, prims, dxinv);
                 , Real dwdx = normal_diff<order>(iv, 0, cls.QW, prims, dxinv);)
#if (AMREX_SPACEDIM >= 2)
    Real dudy = tangent_diff<order>(iv, 0, 1, cls.QU, prims, dxinv);
    Real dvdy = tangent_diff<order>(iv, 0, 1, cls.QV, prims, dxinv);
#endif
#if (AMREX_SPACEDIM == 3)
    Real dudz = tangent_diff<order>(iv, 0, 2, cls.QU, prims, dxinv);
    Real dwdz = tangent_diff<order>(iv, 0, 2, cls.QW, prims, dxinv);
#endif
    Real divu = AMREX_D_TERM(dudx, +dvdy, +dwdz);

    const Real mu = interp<order>(iv, 0, 0, mu_arr);
    const Real xi = interp<order>(iv, 0, 0, xi_arr);
    const Real lam = interp<order>(iv, 0, 0, lam_arr);

    AMREX_D_TERM(Real tauxx = mu * (2.0 * dudx - 2.0 / 3.0 * divu) + xi * divu;
                 , Real tauxy = mu * (dudy + dvdx);
                 , Real tauxz = mu * (dudz + dwdx););

    AMREX_D_TERM(fx(i, j, k, cls.UMX) = -tauxx;, fx(i, j, k, cls.UMY) = -tauxy;
                 , fx(i, j, k, cls.UMZ) = -tauxz;);
    fx(i, j, k, cls.UET) =
        -(AMREX_D_TERM(interp<order>(iv, 0, cls.QU, prims) * tauxx,
                       +interp<order>(iv, 0, cls.QV, prims) * tauxy,
                       +interp<order>(iv, 0, cls.QW, prims) * tauxz)) -
        lam * dTdx;

    // Species transport
    cns_diff_species(iv, 0, prims, rD_arr, dxinv, fx, cls);
  }

#if AMREX_SPACEDIM >= 2
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE void cns_diff_y(
      int i, int j, int k, amrex::Array4<amrex::Real> const& fy,
      amrex::Array4<const amrex::Real> const& prims,
      amrex::Array4<const amrex::Real> const& mu_arr,
      amrex::Array4<const amrex::Real> const& xi_arr,
      amrex::Array4<const amrex::Real> const& lam_arr,
      amrex::Array4<const amrex::Real> const& rD_arr,
      amrex::GpuArray<amrex::Real, amrex::SpaceDim> const& dxinv,
      cls_t const& cls) const noexcept {
    using amrex::Real;
    const amrex::IntVect iv{AMREX_D_DECL(i, j, k)};
    constexpr int order = 2;

    AMREX_D_TERM(, Real dudy = normal_diff<order>(iv, 1, cls.QU, prims, dxinv);
                 Real dvdy = normal_diff<order>(iv, 1, cls.QV, prims, dxinv);
                 , Real dwdy = normal_diff<order>(iv, 1, cls.QW, prims, dxinv);)
    Real dudx = tangent_diff<order>(iv, 1, 0, cls.QU, prims, dxinv);
    Real dvdx = tangent_diff<order>(iv, 1, 0, cls.QV, prims, dxinv);
#if (AMREX_SPACEDIM == 3)
    Real dvdz = tangent_diff<order>(iv, 1, 2, cls.QV, prims, dxinv);
    Real dwdz = tangent_diff<order>(iv, 1, 2, cls.QW, prims, dxinv);
#endif
    Real divu = AMREX_D_TERM(, dudx + dvdy, +dwdz);

    const Real mu = interp<order>(iv, 1, 0, mu_arr);
    const Real xi = interp<order>(iv, 1, 0, xi_arr);

    AMREX_D_TERM(
        , Real tauyy = mu * (2.0 * dvdy - (2.0 / 3.0) * divu) + xi * divu;
        Real tauxy = mu * (dudy + dvdx);, Real tauyz = mu * (dwdy + dvdz);)

    AMREX_D_TERM(, fy(i, j, k, cls.UMX) += -tauxy;
                 fy(i, j, k, cls.UMY) += -tauyy;
                 , fy(i, j, k, cls.UMZ) += -tauyz;);
    Real dTdy = normal_diff<order>(iv, 1, cls.QT, prims, dxinv);
    fy(i, j, k, cls.UET) +=
        -(AMREX_D_TERM(,
                       interp<order>(iv, 0, cls.QU, prims) * tauxy +
                           interp<order>(iv, 0, cls.QV, prims) * tauyy,
                       +interp<order>(iv, 0, cls.QW, prims) * tauyz)) -
        interp<order>(iv, 1, 0, lam_arr) * dTdy;

    // Species transport
  }
#endif
#if AMREX_SPACEDIM == 3
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE void cns_diff_z(
      int i, int j, int k, amrex::Array4<amrex::Real> const& fz,
      amrex::Array4<const amrex::Real> const& prims,
      amrex::Array4<const amrex::Real> const& mu_arr,
      amrex::Array4<const amrex::Real> const& xi_arr,
      amrex::Array4<const amrex::Real> const& lam_arr,
      amrex::Array4<const amrex::Real> const& rD_arr,
      amrex::GpuArray<amrex::Real, amrex::SpaceDim> const& dxinv,
      cls_t const& cls) const noexcept {
    using amrex::Real;
    const amrex::IntVect iv{AMREX_D_DECL(i, j, k)};
    constexpr int order = 2;

    Real dudz = normal_diff<order>(iv, 2, cls.QU, prims, dxinv);
    Real dvdz = normal_diff<order>(iv, 2, cls.QV, prims, dxinv);
    Real dwdz = normal_diff<order>(iv, 2, cls.QW, prims, dxinv);
    Real dudx = tangent_diff<order>(iv, 2, 0, cls.QU, prims, dxinv);
    Real dwdx = tangent_diff<order>(iv, 2, 0, cls.QW, prims, dxinv);
    Real dvdy = tangent_diff<order>(iv, 2, 1, cls.QV, prims, dxinv);
    Real dwdy = tangent_diff<order>(iv, 2, 1, cls.QW, prims, dxinv);
    Real divu = dudx + dvdy + dwdz;

    const Real mu = interp<order>(iv, 2, 0, mu_arr);
    const Real xi = interp<order>(iv, 2, 0, xi_arr);

    Real tauxz = mu * (dudz + dwdx);
    Real tauyz = mu * (dvdz + dwdy);
    Real tauzz = mu * (2.0 * dwdz - (2.0 / 3.0) * divu) + xi * divu;

    fz(i, j, k, cls.UMX) += -tauxz;
    fz(i, j, k, cls.UMY) += -tauyz;
    fz(i, j, k, cls.UMZ) += -tauzz;
    Real dTdz = normal_diff<order>(iv, 2, cls.QT, prims, dxinv);
    fz(i, j, k, cls.UET) += -(interp<order>(iv, 2, cls.QU, prims) * tauxz +
                              interp<order>(iv, 2, cls.QV, prims) * tauyz +
                              interp<order>(iv, 2, cls.QW, prims) * tauzz) -
                            interp<order>(iv, 2, 0, lam_arr) * dTdz;

    // Species transport
  }
#endif
};

#endif
