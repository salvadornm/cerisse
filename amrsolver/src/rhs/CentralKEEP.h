#ifndef CentralKEEP_H_
#define CentralKEEP_H_

#include <Index.h>

///
/// \brief Template class for Central Kinetic Energy and Entropy Preserving
/// scheme
///
/// \param isAD  Artificial dissipation True or False
/// \param isIB  Immersed boundary fluxes
/// \param order Order of accuracy, 2, 4 or 6
///
/// ```
/// {rst}
/// Method implemented from High-order accurate kinetic-energy and entropy
/// preserving (KEEP) schemes on curvilinear grids (10.1016/j.jcp.2021.110482).
/// Divergence split:
/// :math:`\frac{\partial f}{\partial x}\bigg|_i= \frac{\partial}{}`
///
/// Quadratic split:
/// :math:`\frac{\partial fg}{\partial x}\bigg|_i= \frac{\partial}{}`
///
/// Cubic split:
/// :math:`\frac{\partial fgh}{\partial x}\bigg|_i= \frac{\partial}{}`
/// ```
///
template <bool isAD, bool isIB, int order, typename cls_t>
class keep_euler_t {
 public:
  AMREX_GPU_HOST_DEVICE
  keep_euler_t() {
    // initialize coefficients
    coeffs(0, 0) = 1.0;
    coeffs(0, 1) = 0.0;
    coeffs(0, 2) = 0.0;
    coeffs(1, 0) = 4.0 / 3;
    coeffs(1, 1) = -2.0 / 12;
    coeffs(1, 2) = 0.0;
    coeffs(2, 0) = 6.0 / 4;
    coeffs(2, 1) = -6.0 / 20;
    coeffs(2, 2) = 2.0 / 60;
  }

  int order_keep = order;

  // TODO: add copy and move constructors
  // AMREX_GPU_HOST_DEVICE
  // keep_euler_t(const keep_euler_t& other) {}
  // keep_euler_t (keep_euler_t<isAD,isIB,order,cls_t>&& rhs) noexcept {};

  AMREX_GPU_HOST_DEVICE
  ~keep_euler_t() {}

  void inline eflux(const Geometry& geom, const MFIter& mfi,
                    const Array4<Real>& prims, const Array4<Real>& flx,
                    const Array4<Real>& rhs,
                    const cls_t& cls) {
    const GpuArray<Real, AMREX_SPACEDIM> dxinv = geom.InvCellSizeArray();
    const Box& bx  = mfi.growntilebox(0);
    const Box& bxg = mfi.growntilebox(NGHOST);
    const Box& bxgnodal = mfi.grownnodaltilebox(
        -1, 0);  // extent is 0,N_cell+1 in all directions -- -1 means for all
                 // directions. amrex::surroundingNodes(bx) does the same

    int halfsten = order / 2;
    Array1D<Real, 0, 2> order_coeffs{
        coeffs(halfsten - 1, 0), coeffs(halfsten - 1, 1),
        coeffs(halfsten - 1, 2)};  // get coefficients array for current scheme

    // compute interface fluxes at i-1/2, j-1/2, k-1/2
    ParallelFor(bxg, NCONS, [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {rhs(i, j, k, n)=0.0;});

    for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
      ParallelFor(bxgnodal,
                  [=, *this] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    this->flux_dir(i, j, k, dir, order_coeffs, prims, flx, cls);
                  });

      // add flux derivative to rhs = -(fi+1 - fi)/dx = (fi - fi+1)/dx
      ParallelFor(bx, NCONS,
                  [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
                    rhs(i, j, k, n) +=
                        dxinv[dir] * (flx(i, j, k, n) - flx(i + 1, j, k, n));
                  });
    }

    if constexpr (isIB) {
      eflux_ibm();
    }

    if constexpr (isAD) {
      art_dissipation_flux();
    }
  }
  // flux_dir needs to be public only because we need test it. Can improve this.
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE void flux_dir(
      int i, int j, int k, int dir, const Array1D<Real, 0, 2>& coefs,
      const Array4<Real>& prims, const Array4<Real>& flx,
      const cls_t& cls) const {
    for (int n = 0; n < NCONS; n++) {
      flx(i, j, k, n) = 0.0;
    };
    GpuArray<int, 3> vdir = {int(dir == 0), int(dir == 1), int(dir == 2)};

    int i1, j1, i2, j2, k1, k2;
    Real rho1, ux1, uy1, uz1, uu1, ie1, p1;
    Real rho2, ux2, uy2, uz2, uu2, ie2, p2;
    Real massflx, ke;

    // flux at i/j/k-1/2
    for (int l = 1; l <= order / 2; l++) {
      for (int m = 0; m <= l - 1; m++) {
        i1 = i + m * vdir[0];
        j1 = j + m * vdir[1];
        k1 = k + m * vdir[2];
        rho1 = prims(i1, j1, k1, cls.QRHO);
        ux1 = prims(i1, j1, k1, cls.QU);
        uy1 = prims(i1, j1, k1, cls.QV);
        uz1 = prims(i1, j1, k1, cls.QW);
        uu1 = prims(i1, j1, k1, cls.QU + dir);
        ie1 = cls.cv * prims(i1, j1, k1, cls.QT);
        p1 = prims(i1, j1, k1, cls.QPRES);

        i2 = i + (m - l) * vdir[0];
        j2 = j + (m - l) * vdir[1];
        k2 = k + (m - l) * vdir[2];
        rho2 = prims(i2, j2, k2, cls.QRHO);
        ux2 = prims(i2, j2, k2, cls.QU);
        uy2 = prims(i2, j2, k2, cls.QV);
        uz2 = prims(i2, j2, k2, cls.QW);
        uu2 = prims(i2, j2, k2, cls.QU + dir);
        ie2 = cls.cv * prims(i2, j2, k2, cls.QT);
        p2 = prims(i2, j2, k2, cls.QPRES);

        massflx = fgQuad(rho1, rho2, uu1, uu2);
        ke = 0.5_rt * (ux1 * ux2 + uy1 * uy2 + uz1 * uz2);
        flx(i, j, k, cls.URHO) += coefs(l - 1) * massflx;
        flx(i, j, k, cls.UMX) +=
            coefs(l - 1) *
            (fghCubic(rho1, rho2, uu1, uu2, ux1, ux2) + vdir[0] * fDiv(p1, p2));
        flx(i, j, k, cls.UMY) +=
            coefs(l - 1) *
            (fghCubic(rho1, rho2, uu1, uu2, uy1, uy2) + vdir[1] * fDiv(p1, p2));
        flx(i, j, k, cls.UMZ) +=
            coefs(l - 1) *
            (fghCubic(rho1, rho2, uu1, uu2, uz1, uz2) + vdir[2] * fDiv(p1, p2));
        flx(i, j, k, cls.UET) +=
            coefs(l - 1) *
            (massflx * ke + fghCubic(rho1, rho2, uu1, uu2, ie1, ie2) +
             fgDiv(p1, p2, uu1, uu2));
      }
    }
  }

  // immersed boundary ghost point
  void inline eflux_ibm() { amrex::Print() << "IBM eflux" << std::endl; }

  // artificial dissipation flux
  void inline art_dissipation_flux() {
    amrex::Print() << "AD eflux" << std::endl;
  }

  // 2 * standard finite difference coeffs-- need coefficients for all orders
  typedef Array2D<Real, 0, 2, 0, 2> arrCoeff_t;
  arrCoeff_t coeffs;

 private:
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE Real fDiv(Real f, Real fl) const {
    return 0.5_rt * (f + fl);
  }

  // Divergence split
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE Real fgDiv(Real f, Real fl, Real g,
                                                 Real gl) const {
    return 0.5_rt * (f * gl + fl * g);
  }

  // Quadratic split
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE Real fgQuad(Real f, Real fl, Real g,
                                                  Real gl) const {
    return 0.25_rt * (f + fl) * (g + gl);
  }

  // Cubic split
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE Real fghCubic(Real f, Real fl, Real g,
                                                    Real gl, Real h,
                                                    Real hl) const {
    return 0.1250_rt * (f + fl) * (g + gl) * (h + hl);
  }
};

// Immersed boundaries ///////////////////////////////////////////////////
// TODO: the following if results in a compilation error with USE_IBM as all
// branches of constexpr if are evaluated at compile time and IBMultiFab is not
// defined without USE_IBM. This is a problem as we want to use the same code
// for both IBM and non-IBM cases. if constexpr (PROB::do_ibm==1) {
//   IBM::IBMultiFab& mfab = *IBM::ibMFa[level]; // this is ugly
//   // field array
//   auto& ibfab = mfab.get(mfi);
//   auto const& conFabArr  = consmf.array(mfi); // this is a const
//   because.array() returns a const but we can still modify conFab as consmf
//   input argument is not const auto const& primFabArr = primsmf.array(mfi);
//   auto const& ibFabArr = mfab.array(mfi);

//   IBM::IBMultiFab& mfab = *ibMFa[level];
//   IBM::ib.computeGPs(consfab, primsfab, ibfab, cls);
// }

#endif