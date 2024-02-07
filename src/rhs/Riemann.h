#ifndef Riemann_H_
#define Riemann_H_

#include <AMReX_FArrayBox.H>
#include <CNS.h>

///
/// \brief Template class for Riemann solvers
///
/// \param iOption Flux vector splitting (0 - Global lax-friedrichs)
///
/// ```
/// {rst}
///
/// :math:`f_{i+1/2}= `
/// ```
///
template <bool iOption, typename cls_t>
class riemann_t {
 public:
  AMREX_GPU_HOST_DEVICE
  riemann_t() {}

  AMREX_GPU_HOST_DEVICE
  ~riemann_t() {}

  void inline eflux(const Geometry& geom, const MFIter& mfi,
                    const Array4<Real>& prims, const Array4<Real>& flx,
                    const Array4<Real>& rhs, const cls_t& cls) {
    const GpuArray<Real, AMREX_SPACEDIM> dxinv = geom.InvCellSizeArray();
    const Box& bx = mfi.tilebox();
    const Box& bxg = mfi.growntilebox(cls_t::NGHOST);
    const Box& bxn = mfi.grownnodaltilebox(-1, 0);  // 0,N+1 all directions

    // zero rhs
    ParallelFor(bxg, cls_t::NCONS,
                [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
                  rhs(i, j, k, n) = 0.0;
                });

    const Box& bxg1 = amrex::grow(bx, 1);
    FArrayBox slopef(bxg1, cls_t::NCONS, The_Async_Arena());
    const Array4<Real>& slope = slopef.array();

    // x-direction
    int cdir = 0;
    const Box& xslpbx = amrex::grow(bx, cdir, 1);
    amrex::ParallelFor(
        bx, [=, *this] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
          this->cns_slope_x(i, j, k, slope, prims, cls);
        });
    const Box& xflxbx = amrex::surroundingNodes(bx, cdir);
    amrex::ParallelFor(
        xflxbx, [=, *this] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
          this->cns_riemann_x(i, j, k, flx, slope, prims, cls);
        });
    // add x flux derivative to rhs = -(fi+1 - fi)/dx = (fi - fi+1)/dx
    ParallelFor(bx, cls_t::NCONS,
                [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
                  rhs(i, j, k, n) +=
                      dxinv[cdir] * (flx(i, j, k, n) - flx(i + 1, j, k, n));
                });

    // y-direction
    cdir = 1;
    const Box& yslpbx = amrex::grow(bx, cdir, 1);
    amrex::ParallelFor(
        yslpbx, [=, *this] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
          this->cns_slope_y(i, j, k, slope, prims, cls);
        });
    const Box& yflxbx = amrex::surroundingNodes(bx, cdir);
    amrex::ParallelFor(
        yflxbx, [=, *this] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
          for (int n = 0; n < cls_t::NCONS; n++) {
            flx(i, j, k, n) = 0.0;
          };
          this->cns_riemann_y(i, j, k, flx, slope, prims, cls);
        });
    // add y flux derivative to rhs = -(fi+1 - fi)/dy = (fi - fi+1)/dy
    ParallelFor(bx, cls_t::NCONS,
                [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
                  rhs(i, j, k, n) +=
                      dxinv[cdir] * (flx(i, j, k, n) - flx(i, j + 1, k, n));
                });

    // z-direction
    cdir = 2;
    const Box& zslpbx = amrex::grow(bx, cdir, 1);
    amrex::ParallelFor(
        zslpbx, [=, *this] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
          this->cns_slope_z(i, j, k, slope, prims, cls);
        });
    const Box& zflxbx = amrex::surroundingNodes(bx, cdir);
    amrex::ParallelFor(
        zflxbx, [=, *this] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
          for (int n = 0; n < cls_t::NCONS; n++) {
            flx(i, j, k, n) = 0.0;
          };
          this->cns_riemann_z(i, j, k, flx, slope, prims, cls);
        });
    // add z flux derivative to rhs = -(fi+1 - fi)/dz = (fi - fi+1)/dz
    ParallelFor(bx, cls_t::NCONS,
                [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
                  rhs(i, j, k, n) +=
                      dxinv[cdir] * (flx(i, j, k, n) - flx(i, j, k + 1, n));
                });
  };

  AMREX_GPU_DEVICE AMREX_FORCE_INLINE Real limiter(Real dlft, Real drgt) const {
    Real dcen = Real(0.5) * (dlft + drgt);
    Real dsgn = Math::copysign(Real(1.0), dcen);
    Real slop = Real(2.0) * min(Math::abs(dlft), Math::abs(drgt));
    Real dlim = (dlft * drgt >= Real(0.0)) ? slop : Real(0.0);
    return dsgn * min(dlim, Math::abs(dcen));
  }

  AMREX_GPU_DEVICE AMREX_FORCE_INLINE void riemann_prob(
      const Real gamma, const Real smallp, const Real /*smallr*/, const Real rl,
      const Real ul, const Real pl, const Real ut1l, const Real ut2l,
      const Real rr, const Real ur, const Real pr, const Real ut1r,
      const Real ut2r, Real& flxrho, Real& flxu, Real& flxut, Real& flxutt,
      Real& flxe) const {
    constexpr Real weakwv = Real(1.e-3);
    constexpr Real small = Real(1.e-6);

    Real clsql = gamma * pl * rl;
    Real clsqr = gamma * pr * rr;
    Real wl = std::sqrt(clsql);
    Real wr = std::sqrt(clsqr);
    Real cleft = wl / rl;
    Real cright = wr / rr;
    Real ccsmall = small * (cleft + cright);

    Real pstar = (wl * pr + wr * pl - wr * wl * (ur - ul)) / (wl + wr);
    pstar = max(pstar, smallp);
    Real pstnm1 = pstar;

    Real wlsq = (Real(0.5) * (gamma - Real(1.)) * (pstar + pl) + pstar) * rl;
    Real wrsq = (Real(0.5) * (gamma - Real(1.)) * (pstar + pr) + pstar) * rr;

    wl = std::sqrt(wlsq);
    wr = std::sqrt(wrsq);
    Real ustarp = ul - (pstar - pl) / wl;
    Real ustarm = ur + (pstar - pr) / wr;

    pstar = (wl * pr + wr * pl - wr * wl * (ur - ul)) / (wl + wr);
    pstar = max(pstar, smallp);

    Real ustar;
    for (int iter = 0; iter < 3; ++iter) {
      wlsq = (Real(0.5) * (gamma - Real(1.)) * (pstar + pl) + pstar) * rl;
      wrsq = (Real(0.5) * (gamma - Real(1.)) * (pstar + pr) + pstar) * rr;

      wl = Real(1.) / std::sqrt(wlsq);
      wr = Real(1.) / std::sqrt(wrsq);

      Real ustnm1 = ustarm;
      Real ustnp1 = ustarp;

      ustarm = ur - (pr - pstar) * wr;
      ustarp = ul + (pl - pstar) * wl;

      Real dpditer = Math::abs(pstnm1 - pstar);
      Real zp = Math::abs(ustarp - ustnp1);
      if (zp - weakwv * cleft < Real(0.0)) {
        zp = dpditer * wl;
      }
      Real zm = Math::abs(ustarm - ustnm1);
      if (zm - weakwv * cright < Real(0.0)) {
        zm = dpditer * wr;
      }

      Real zz = zp + zm;
      Real denom = dpditer / max(zz, ccsmall);
      pstnm1 = pstar;
      pstar = pstar - denom * (ustarm - ustarp);
      pstar = max(pstar, smallp);
      ustar = Real(0.5) * (ustarm + ustarp);
    }

    Real ro, uo, po, sgnm, utrans1, utrans2;
    if (ustar > Real(0.)) {
      ro = rl;
      uo = ul;
      po = pl;
      sgnm = Real(1.);
      utrans1 = ut1l;
      utrans2 = ut2l;
    } else if (ustar < Real(0.)) {
      ro = rr;
      uo = ur;
      po = pr;
      sgnm = Real(-1.);
      utrans1 = ut1r;
      utrans2 = ut2r;
    } else {
      uo = Real(0.5) * (ur + ul);
      po = Real(0.5) * (pr + pl);
      ro = Real(2.) * (rl * rr) / (rl + rr);
      sgnm = Real(1.);
      utrans1 = Real(0.5) * (ut1l + ut1r);
      utrans2 = Real(0.5) * (ut2l + ut2r);
    }
    Real wosq = (Real(0.5) * (gamma - Real(1.)) * (pstar + po) + pstar) * ro;
    Real co = std::sqrt(gamma * po / ro);
    Real wo = std::sqrt(wosq);
    Real dpjmp = pstar - po;
    Real rstar = ro / (Real(1.) - ro * dpjmp / wosq);
    Real cstar = std::sqrt(gamma * pstar / rstar);
    Real spout = co - sgnm * uo;
    Real spin = cstar - sgnm * uo;
    if (pstar >= po) {
      spin = wo / ro - sgnm * uo;
      spout = spin;
    }
    Real ss = max(spout - spin, spout + spin);
    Real frac = Real(0.5) * (Real(1.) + (spin + spout) / max(ss, ccsmall));

    Real rgdnv, ugdnv, pgdnv;
    if (spout < Real(0.)) {
      rgdnv = ro;
      ugdnv = uo;
      pgdnv = po;
    } else if (spin >= Real(0.)) {
      rgdnv = rstar;
      ugdnv = ustar;
      pgdnv = pstar;
    } else {
      rgdnv = frac * rstar + (Real(1.) - frac) * ro;
      ugdnv = frac * ustar + (Real(1.) - frac) * uo;
      pgdnv = frac * pstar + (Real(1.) - frac) * po;
    }

    flxrho = rgdnv * ugdnv;
    flxu = rgdnv * ugdnv * ugdnv + pgdnv;
    flxut = rgdnv * ugdnv * utrans1;
    flxutt = rgdnv * ugdnv * utrans2;
    flxe =
        ugdnv * (Real(0.5) * rgdnv *
                     (ugdnv * ugdnv + utrans1 * utrans1 + utrans2 * utrans2) +
                 pgdnv / (gamma - Real(1.)) + pgdnv);
  }

  AMREX_GPU_DEVICE AMREX_FORCE_INLINE void cns_slope_x(
      int i, int j, int k, const Array4<Real>& dq, const Array4<Real>& q,
      const cls_t& cls) const {
    Real cspeed = sqrt(cls.gamma * cls.Rspec * q(i, j, k, cls.QT)) + 2.e-40;
    Real dlft = Real(0.5) *
                    (q(i, j, k, cls.QPRES) - q(i - 1, j, k, cls.QPRES)) /
                    cspeed -
                Real(0.5) * q(i, j, k, cls.QRHO) *
                    (q(i, j, k, cls.QU) - q(i - 1, j, k, cls.QU));
    Real drgt = Real(0.5) *
                    (q(i + 1, j, k, cls.QPRES) - q(i, j, k, cls.QPRES)) /
                    cspeed -
                Real(0.5) * q(i, j, k, cls.QRHO) *
                    (q(i + 1, j, k, cls.QU) - q(i, j, k, cls.QU));
    Real d0 = limiter(dlft, drgt);

    Real cs2 = cspeed * cspeed;
    dlft = (q(i, j, k, cls.QRHO) - q(i - 1, j, k, cls.QRHO)) -
           (q(i, j, k, cls.QPRES) - q(i - 1, j, k, cls.QPRES)) / cs2;
    drgt = (q(i + 1, j, k, cls.QRHO) - q(i, j, k, cls.QRHO)) -
           (q(i + 1, j, k, cls.QPRES) - q(i, j, k, cls.QPRES)) / cs2;
    Real d1 = limiter(dlft, drgt);

    dlft = Real(0.5) * (q(i, j, k, cls.QPRES) - q(i - 1, j, k, cls.QPRES)) /
               cspeed +
           Real(0.5) * q(i, j, k, cls.QRHO) *
               (q(i, j, k, cls.QU) - q(i - 1, j, k, cls.QU));
    drgt = Real(0.5) * (q(i + 1, j, k, cls.QPRES) - q(i, j, k, cls.QPRES)) /
               cspeed +
           Real(0.5) * q(i, j, k, cls.QRHO) *
               (q(i + 1, j, k, cls.QU) - q(i, j, k, cls.QU));
    Real d2 = limiter(dlft, drgt);

    dlft = q(i, j, k, cls.QV) - q(i - 1, j, k, cls.QV);
    drgt = q(i + 1, j, k, cls.QV) - q(i, j, k, cls.QV);
    Real d3 = limiter(dlft, drgt);

    dlft = q(i, j, k, cls.QW) - q(i - 1, j, k, cls.QW);
    drgt = q(i + 1, j, k, cls.QW) - q(i, j, k, cls.QW);
    Real d4 = limiter(dlft, drgt);

    dq(i, j, k, 0) = d0;
    dq(i, j, k, 1) = d1;
    dq(i, j, k, 2) = d2;
    dq(i, j, k, 3) = d3;
    dq(i, j, k, 4) = d4;
  }

  AMREX_GPU_DEVICE AMREX_FORCE_INLINE void cns_slope_y(
      int i, int j, int k, Array4<Real> const& dq, Array4<Real const> const& q,
      cls_t const& cls) const {
    Real cspeed = sqrt(cls.gamma * cls.Rspec * q(i, j, k, cls.QT)) + 1.e-40;
    Real dlft = Real(0.5) *
                    (q(i, j, k, cls.QPRES) - q(i, j - 1, k, cls.QPRES)) /
                    cspeed -
                Real(0.5) * q(i, j, k, cls.QRHO) *
                    (q(i, j, k, cls.QV) - q(i, j - 1, k, cls.QV));
    Real drgt = Real(0.5) *
                    (q(i, j + 1, k, cls.QPRES) - q(i, j, k, cls.QPRES)) /
                    cspeed -
                Real(0.5) * q(i, j, k, cls.QRHO) *
                    (q(i, j + 1, k, cls.QV) - q(i, j, k, cls.QV));
    Real d0 = limiter(dlft, drgt);

    Real cs2 = cspeed * cspeed;
    dlft = (q(i, j, k, cls.QRHO) - q(i, j - 1, k, cls.QRHO)) -
           (q(i, j, k, cls.QPRES) - q(i, j - 1, k, cls.QPRES)) / cs2;
    drgt = (q(i, j + 1, k, cls.QRHO) - q(i, j, k, cls.QRHO)) -
           (q(i, j + 1, k, cls.QPRES) - q(i, j, k, cls.QPRES)) / cs2;
    Real d1 = limiter(dlft, drgt);

    dlft = Real(0.5) * (q(i, j, k, cls.QPRES) - q(i, j - 1, k, cls.QPRES)) /
               cspeed +
           Real(0.5) * q(i, j, k, cls.QRHO) *
               (q(i, j, k, cls.QV) - q(i, j - 1, k, cls.QV));
    drgt = Real(0.5) * (q(i, j + 1, k, cls.QPRES) - q(i, j, k, cls.QPRES)) /
               cspeed +
           Real(0.5) * q(i, j, k, cls.QRHO) *
               (q(i, j + 1, k, cls.QV) - q(i, j, k, cls.QV));
    Real d2 = limiter(dlft, drgt);

    dlft = q(i, j, k, cls.QU) - q(i, j - 1, k, cls.QU);
    drgt = q(i, j + 1, k, cls.QU) - q(i, j, k, cls.QU);
    Real d3 = limiter(dlft, drgt);

    dlft = q(i, j, k, cls.QW) - q(i, j - 1, k, cls.QW);
    drgt = q(i, j + 1, k, cls.QW) - q(i, j, k, cls.QW);
    Real d4 = limiter(dlft, drgt);

    dq(i, j, k, 0) = d0;
    dq(i, j, k, 1) = d1;
    dq(i, j, k, 2) = d2;
    dq(i, j, k, 3) = d3;
    dq(i, j, k, 4) = d4;
  }

  AMREX_GPU_DEVICE AMREX_FORCE_INLINE void cns_slope_z(
      int i, int j, int k, Array4<Real> const& dq, Array4<Real const> const& q,
      cls_t const& cls) const {
    Real cspeed = sqrt(cls.gamma * cls.Rspec * q(i, j, k, cls.QT)) + 1.e-40;
    Real dlft = Real(0.5) *
                    (q(i, j, k, cls.QPRES) - q(i, j, k - 1, cls.QPRES)) /
                    cspeed -
                Real(0.5) * q(i, j, k, cls.QRHO) *
                    (q(i, j, k, cls.QW) - q(i, j, k - 1, cls.QW));
    Real drgt = Real(0.5) *
                    (q(i, j, k + 1, cls.QPRES) - q(i, j, k, cls.QPRES)) /
                    cspeed -
                Real(0.5) * q(i, j, k, cls.QRHO) *
                    (q(i, j, k + 1, cls.QW) - q(i, j, k, cls.QW));
    Real d0 = limiter(dlft, drgt);

    Real cs2 = cspeed * cspeed;
    dlft = (q(i, j, k, cls.QRHO) - q(i, j, k - 1, cls.QRHO)) -
           (q(i, j, k, cls.QPRES) - q(i, j, k - 1, cls.QPRES)) / cs2;
    drgt = (q(i, j, k + 1, cls.QRHO) - q(i, j, k, cls.QRHO)) -
           (q(i, j, k + 1, cls.QPRES) - q(i, j, k, cls.QPRES)) / cs2;
    Real d1 = limiter(dlft, drgt);

    dlft = Real(0.5) * (q(i, j, k, cls.QPRES) - q(i, j, k - 1, cls.QPRES)) /
               cspeed +
           Real(0.5) * q(i, j, k, cls.QRHO) *
               (q(i, j, k, cls.QW) - q(i, j, k - 1, cls.QW));
    drgt = Real(0.5) * (q(i, j, k + 1, cls.QPRES) - q(i, j, k, cls.QPRES)) /
               cspeed +
           Real(0.5) * q(i, j, k, cls.QRHO) *
               (q(i, j, k + 1, cls.QW) - q(i, j, k, cls.QW));
    Real d2 = limiter(dlft, drgt);

    dlft = q(i, j, k, cls.QU) - q(i, j, k - 1, cls.QU);
    drgt = q(i, j, k + 1, cls.QU) - q(i, j, k, cls.QU);
    Real d3 = limiter(dlft, drgt);

    dlft = q(i, j, k, cls.QV) - q(i, j, k - 1, cls.QV);
    drgt = q(i, j, k + 1, cls.QV) - q(i, j, k, cls.QV);
    Real d4 = limiter(dlft, drgt);

    dq(i, j, k, 0) = d0;
    dq(i, j, k, 1) = d1;
    dq(i, j, k, 2) = d2;
    dq(i, j, k, 3) = d3;
    dq(i, j, k, 4) = d4;
  }

  AMREX_GPU_DEVICE AMREX_FORCE_INLINE void cns_riemann_x(
      int i, int j, int k, Array4<Real> const& fx, Array4<Real> const& dq,
      Array4<Real> const& q, cls_t const& cls) const {
    Real cspeed = sqrt(cls.gamma * cls.Rspec * q(i - 1, j, k, cls.QT)) + 1.e-40;
    Real rl = q(i - 1, j, k, cls.QRHO) +
              Real(0.5) * ((dq(i - 1, j, k, 0) + dq(i - 1, j, k, 2)) / cspeed +
                           dq(i - 1, j, k, 1));
    rl = max(rl, 1.0e-40);
    Real ul = q(i - 1, j, k, cls.QU) +
              Real(0.5) * ((dq(i - 1, j, k, 2) - dq(i - 1, j, k, 0)) /
                           q(i - 1, j, k, cls.QRHO));
    Real pl = q(i - 1, j, k, cls.QPRES) +
              Real(0.5) * (dq(i - 1, j, k, 0) + dq(i - 1, j, k, 2)) * cspeed;
    pl = max(pl, 1.0e-40);
    Real ut1l = q(i - 1, j, k, cls.QV) + Real(0.5) * dq(i - 1, j, k, 3);
    Real ut2l = q(i - 1, j, k, cls.QW) + Real(0.5) * dq(i - 1, j, k, 4);

    cspeed = sqrt(cls.gamma * cls.Rspec * q(i, j, k, cls.QT)) + 1.e-40;
    Real rr = q(i, j, k, cls.QRHO) -
              Real(0.5) *
                  ((dq(i, j, k, 0) + dq(i, j, k, 2)) / cspeed + dq(i, j, k, 1));
    rr = max(rr, 1.0e-40);
    Real ur =
        q(i, j, k, cls.QU) -
        Real(0.5) * ((dq(i, j, k, 2) - dq(i, j, k, 0)) / q(i, j, k, cls.QRHO));
    Real pr = q(i, j, k, cls.QPRES) -
              Real(0.5) * (dq(i, j, k, 0) + dq(i, j, k, 2)) * cspeed;
    pr = max(pr, 1.0e-40);
    Real ut1r = q(i, j, k, cls.QV) - Real(0.5) * dq(i, j, k, 3);
    Real ut2r = q(i, j, k, cls.QW) - Real(0.5) * dq(i, j, k, 4);

    riemann_prob(cls.gamma, 1.0e-40, 1.0e-40, rl, ul, pl, ut1l, ut2l, rr, ur,
                 pr, ut1r, ut2r, fx(i, j, k, cls.URHO), fx(i, j, k, cls.UMX),
                 fx(i, j, k, cls.UMY), fx(i, j, k, cls.UMZ),
                 fx(i, j, k, cls.UET));
  }

  AMREX_GPU_DEVICE AMREX_FORCE_INLINE void cns_riemann_y(
      int i, int j, int k, Array4<Real> const& fy, Array4<Real const> const& dq,
      Array4<Real const> const& q, cls_t const& cls) const {
    Real cspeed = sqrt(cls.gamma * cls.Rspec * q(i, j - 1, k, cls.QT)) + 1.e-40;
    Real rl = q(i, j - 1, k, cls.QRHO) +
              Real(0.5) * ((dq(i, j - 1, k, 0) + dq(i, j - 1, k, 2)) / cspeed +
                           dq(i, j - 1, k, 1));
    rl = max(rl, 1.0e-40);
    Real ul = q(i, j - 1, k, cls.QV) +
              Real(0.5) * ((dq(i, j - 1, k, 2) - dq(i, j - 1, k, 0)) /
                           q(i, j - 1, k, cls.QRHO));
    Real pl = q(i, j - 1, k, cls.QPRES) +
              Real(0.5) * (dq(i, j - 1, k, 0) + dq(i, j - 1, k, 2)) * cspeed;
    pl = max(pl, 1.0e-40);
    Real ut1l = q(i, j - 1, k, cls.QU) + Real(0.5) * dq(i, j - 1, k, 3);
    Real ut2l = q(i, j - 1, k, cls.QW) + Real(0.5) * dq(i, j - 1, k, 4);

    cspeed = sqrt(cls.gamma * cls.Rspec * q(i, j, k, cls.QT)) + 1.e-40;
    Real rr = q(i, j, k, cls.QRHO) -
              Real(0.5) *
                  ((dq(i, j, k, 0) + dq(i, j, k, 2)) / cspeed + dq(i, j, k, 1));
    rr = max(rr, 1.0e-40);
    Real ur =
        q(i, j, k, cls.QV) -
        Real(0.5) * ((dq(i, j, k, 2) - dq(i, j, k, 0)) / q(i, j, k, cls.QRHO));
    Real pr = q(i, j, k, cls.QPRES) -
              Real(0.5) * (dq(i, j, k, 0) + dq(i, j, k, 2)) * cspeed;
    pr = max(pr, 1.0e-40);
    Real ut1r = q(i, j, k, cls.QU) - Real(0.5) * dq(i, j, k, 3);
    Real ut2r = q(i, j, k, cls.QW) - Real(0.5) * dq(i, j, k, 4);

    riemann_prob(cls.gamma, 1.0e-40, 1.0e-40, rl, ul, pl, ut1l, ut2l, rr, ur,
                 pr, ut1r, ut2r, fy(i, j, k, cls.URHO), fy(i, j, k, cls.UMY),
                 fy(i, j, k, cls.UMX), fy(i, j, k, cls.UMZ),
                 fy(i, j, k, cls.UET));
  }

  AMREX_GPU_DEVICE AMREX_FORCE_INLINE void cns_riemann_z(
      int i, int j, int k, Array4<Real> const& fz, Array4<Real const> const& dq,
      Array4<Real const> const& q, cls_t const& cls) const {
    Real cspeed = sqrt(cls.gamma * cls.Rspec * q(i, j, k + 1, cls.QT)) + 1.e-40;
    Real rl = q(i, j, k - 1, cls.QRHO) +
              Real(0.5) * ((dq(i, j, k - 1, 0) + dq(i, j, k - 1, 2)) / cspeed +
                           dq(i, j, k - 1, 1));
    rl = max(rl, 1.0e-40);
    Real ul = q(i, j, k - 1, cls.QW) +
              Real(0.5) * ((dq(i, j, k - 1, 2) - dq(i, j, k - 1, 0)) /
                           q(i, j, k - 1, cls.QRHO));
    Real pl = q(i, j, k - 1, cls.QPRES) +
              Real(0.5) * (dq(i, j, k - 1, 0) + dq(i, j, k - 1, 2)) * cspeed;
    pl = max(pl, 1.0e-40);
    Real ut1l = q(i, j, k - 1, cls.QU) + Real(0.5) * dq(i, j, k - 1, 3);
    Real ut2l = q(i, j, k - 1, cls.QV) + Real(0.5) * dq(i, j, k - 1, 4);

    cspeed = sqrt(cls.gamma * cls.Rspec * q(i, j, k, cls.QT)) + 1.e-40;
    Real rr = q(i, j, k, cls.QRHO) -
              Real(0.5) *
                  ((dq(i, j, k, 0) + dq(i, j, k, 2)) / cspeed + dq(i, j, k, 1));
    rr = max(rr, 1.0e-40);
    Real ur =
        q(i, j, k, cls.QW) -
        Real(0.5) * ((dq(i, j, k, 2) - dq(i, j, k, 0)) / q(i, j, k, cls.QRHO));
    Real pr = q(i, j, k, cls.QPRES) -
              Real(0.5) * (dq(i, j, k, 0) + dq(i, j, k, 2)) * cspeed;
    pr = max(pr, 1.0e-40);
    Real ut1r = q(i, j, k, cls.QU) - Real(0.5) * dq(i, j, k, 3);
    Real ut2r = q(i, j, k, cls.QV) - Real(0.5) * dq(i, j, k, 4);

    riemann_prob(cls.gamma, 1.0e-40, 1.0e-40, rl, ul, pl, ut1l, ut2l, rr, ur,
                 pr, ut1r, ut2r, fz(i, j, k, cls.URHO), fz(i, j, k, cls.UMZ),
                 fz(i, j, k, cls.UMX), fz(i, j, k, cls.UMY),
                 fz(i, j, k, cls.UET));
  }
};

#endif