#ifndef CENTRAL_KEEP_SPEC_H_
#define CENTRAL_KEEP_SPEC_H_

#include <AMReX_Array.H>
#include <AMReX_FArrayBox.H>

/// Multispecies extension of the central KEEP flux.
///
/// CentralKEEP.h assumes a calorically perfect single-species gas: it writes a
/// standalone density flux at URHO, uses cv*T for internal energy, and does not
/// write rho*Y fluxes. In a PelePhysics closure URHO aliases UFS, so that flux
/// is incomplete and inconsistent. This variant instead:
///
///  * obtains e(T,Y,rho) from the selected EOS;
///  * transports every rho*Y_k with the same split mass flux; and
///  * makes sum_k F(rho*Y_k) equal the split mass flux when sum_k Y_k = 1.
///
/// The momentum and total-energy split forms otherwise match CentralKEEP.h.
/// The arithmetic face mass fractions provide a conservative,
/// composition-aware baseline intended to reduce pressure-equilibrium errors;
/// a formally
/// entropy-conservative multicomponent average can replace species_average()
/// later without changing the solver interface.
template <bool isAD, bool isIB, int order, typename cls_t>
class keep_euler_spec_t {
 public:
  AMREX_GPU_HOST_DEVICE
  keep_euler_spec_t() {
    coeffs(0, 0) = 1.0;
    coeffs(0, 1) = 0.0;
    coeffs(0, 2) = 0.0;
    coeffs(1, 0) = 4.0 / 3.0;
    coeffs(1, 1) = -2.0 / 12.0;
    coeffs(1, 2) = 0.0;
    coeffs(2, 0) = 6.0 / 4.0;
    coeffs(2, 1) = -6.0 / 20.0;
    coeffs(2, 2) = 2.0 / 60.0;
  }

  AMREX_GPU_HOST_DEVICE ~keep_euler_spec_t() = default;

  int order_keep = order;

  void inline eflux(
      const Geometry& /*geom*/, const MFIter& mfi,
      const Array4<Real>& prims,
      std::array<FArrayBox*, AMREX_SPACEDIM> const& flxt,
      const Array4<Real>& /*cons*/, const cls_t* cls) {
    static_assert(NUM_SPECIES > 1,
                  "keep_euler_spec_t requires a multispecies closure");
    static_assert(order == 2 || order == 4 || order == 6,
                  "KEEP order must be 2, 4, or 6");

    const Box& bxgnodal = mfi.grownnodaltilebox(-1, 0);
    const int halfsten = order / 2;
    const Array1D<Real, 0, 2> order_coeffs{
        coeffs(halfsten - 1, 0), coeffs(halfsten - 1, 1),
        coeffs(halfsten - 1, 2)};

    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
      auto const& flx = flxt[dir]->array();
      ParallelFor(bxgnodal,
                  [=, *this] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    flux_dir(i, j, k, dir, order_coeffs, prims, flx, cls);
                  });
    }

    if constexpr (isIB) { eflux_ibm(); }
    if constexpr (isAD) { art_dissipation_flux(); }
  }

  AMREX_GPU_DEVICE AMREX_FORCE_INLINE
  void flux_dir(int i, int j, int k, int dir,
                const Array1D<Real, 0, 2>& coefs,
                const Array4<Real>& q, const Array4<Real>& flx,
                const cls_t* cls) const {
    for (int n = 0; n < cls_t::NCONS; ++n) { flx(i, j, k, n) = Real(0.0); }

    const GpuArray<int, 3> vdir = {
        int(dir == 0), int(dir == 1), int(dir == 2)};

    for (int l = 1; l <= order / 2; ++l) {
      for (int m = 0; m <= l - 1; ++m) {
        const int i1 = i + m*vdir[0];
        const int j1 = j + m*vdir[1];
        const int k1 = k + m*vdir[2];
        const int i2 = i + (m-l)*vdir[0];
        const int j2 = j + (m-l)*vdir[1];
        const int k2 = k + (m-l)*vdir[2];

        const Real rho1 = q(i1,j1,k1,cls_t::QRHO);
        const Real rho2 = q(i2,j2,k2,cls_t::QRHO);
        const Real un1 = q(i1,j1,k1,cls_t::QU + dir);
        const Real un2 = q(i2,j2,k2,cls_t::QU + dir);
        const Real ux1 = q(i1,j1,k1,cls_t::QU);
        const Real ux2 = q(i2,j2,k2,cls_t::QU);
        const Real uy1 = q(i1,j1,k1,cls_t::QV);
        const Real uy2 = q(i2,j2,k2,cls_t::QV);
        const Real uz1 = q(i1,j1,k1,cls_t::QW);
        const Real uz2 = q(i2,j2,k2,cls_t::QW);
        const Real p1 = q(i1,j1,k1,cls_t::QPRES);
        const Real p2 = q(i2,j2,k2,cls_t::QPRES);

        Real Y1[NUM_SPECIES];
        Real Y2[NUM_SPECIES];
        for (int ns = 0; ns < NUM_SPECIES; ++ns) {
          Y1[ns] = q(i1,j1,k1,cls_t::QFS + ns);
          Y2[ns] = q(i2,j2,k2,cls_t::QFS + ns);
        }

        Real e1 = Real(0.0);
        Real e2 = Real(0.0);
        cls->TYR2E(q(i1,j1,k1,cls_t::QT), Y1, rho1, e1);
        cls->TYR2E(q(i2,j2,k2,cls_t::QT), Y2, rho2, e2);

        const Real mass_flux = fg_quad(rho1, rho2, un1, un2);
        const Real kinetic_cross = Real(0.5)*
            (ux1*ux2 + uy1*uy2 + uz1*uz2);
        const Real weight = coefs(l-1);

        flx(i,j,k,cls_t::UMX) += weight*
            (fgh_cubic(rho1,rho2,un1,un2,ux1,ux2) +
             Real(vdir[0])*f_div(p1,p2));
        flx(i,j,k,cls_t::UMY) += weight*
            (fgh_cubic(rho1,rho2,un1,un2,uy1,uy2) +
             Real(vdir[1])*f_div(p1,p2));
        flx(i,j,k,cls_t::UMZ) += weight*
            (fgh_cubic(rho1,rho2,un1,un2,uz1,uz2) +
             Real(vdir[2])*f_div(p1,p2));
        flx(i,j,k,cls_t::UET) += weight*
            (mass_flux*kinetic_cross +
             fgh_cubic(rho1,rho2,un1,un2,e1,e2) +
             fg_div(p1,p2,un1,un2));

        for (int ns = 0; ns < NUM_SPECIES; ++ns) {
          flx(i,j,k,cls_t::UFS + ns) +=
              weight*mass_flux*species_average(Y1[ns],Y2[ns]);
        }
      }
    }
  }

  void inline eflux_ibm() {
    amrex::Abort("keep_euler_spec_t: immersed-boundary flux is not implemented");
  }

  void inline art_dissipation_flux() {
    amrex::Abort("keep_euler_spec_t: artificial dissipation is not implemented");
  }

 private:
  using coeff_array_t = Array2D<Real, 0, 2, 0, 2>;
  coeff_array_t coeffs;

  AMREX_GPU_DEVICE AMREX_FORCE_INLINE
  Real f_div(Real a, Real b) const { return Real(0.5)*(a+b); }

  AMREX_GPU_DEVICE AMREX_FORCE_INLINE
  Real fg_div(Real a, Real b, Real c, Real d) const {
    return Real(0.5)*(a*d+b*c);
  }

  AMREX_GPU_DEVICE AMREX_FORCE_INLINE
  Real fg_quad(Real a, Real b, Real c, Real d) const {
    return Real(0.25)*(a+b)*(c+d);
  }

  AMREX_GPU_DEVICE AMREX_FORCE_INLINE
  Real fgh_cubic(Real a, Real b, Real c, Real d, Real e, Real f) const {
    return Real(0.125)*(a+b)*(c+d)*(e+f);
  }

  AMREX_GPU_DEVICE AMREX_FORCE_INLINE
  Real species_average(Real yl, Real yr) const { return Real(0.5)*(yl+yr); }
};

#endif
