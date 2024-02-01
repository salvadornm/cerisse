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
template <bool iOption, typename closures>
class riemann_t {
  public: 
  AMREX_GPU_HOST_DEVICE
  riemann_t() {}
  ~riemann_t() {}

  void inline eflux(const Geometry& geom, const MFIter& mfi,
                    const Array4<Real>& prims, const Array4<Real>& ivars,
                    const Array4<Real>& rhs,
                    const closures& cls) {
    const GpuArray<Real, AMREX_SPACEDIM> dxinv = geom.InvCellSizeArray();
    const Box& bx = mfi.tilebox();
    const Box& bxg = mfi.growntilebox(NGHOST);
    const Box& bxn = mfi.grownnodaltilebox(-1, 0);  // 0,N+1 all directions

    // zero rhs
    ParallelFor(bxg, NCONS, [=] AMREX_GPU_DEVICE(int i, int j, int k, int
    n) noexcept {rhs(i, j, k, n)=0.0;});

    AMREX_D_TERM(auto const& nfabfx = numflxmf[0].array(mfi);
                 , auto const& nfabfy = numflxmf[1].array(mfi);
                 , auto const& nfabfz = numflxmf[2].array(mfi););

    const Box& bxg1 = amrex::grow(bx, 1);
    slopetmp.resize(bxg1, NCONS);
    Elixir slopeeli = slopetmp.elixir();

    FArrayBox slopef(bxg, NCONS, The_Async_Arena());
    auto const& slope = slope.array();

    // x-direction
    int cdir = 0;
    const Box& xslpbx = amrex::grow(bx, cdir, 1);
    amrex::ParallelFor(xslpbx,
                       [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                         cns_slope_x(i, j, k, slope, prims, cls);
                       });
    const Box& xflxbx = amrex::surroundingNodes(bx, cdir);
    amrex::ParallelFor(xflxbx,
                       [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                         cns_riemann_x(i, j, k, nfabfx, slope, prims, cls);
                       });
    
    

    // y-direction
    cdir = 1;
    const Box& yslpbx = amrex::grow(bx, cdir, 1);
    amrex::ParallelFor(yslpbx,
                       [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                         cns_slope_y(i, j, k, slope, prims, cls);
                       });
    const Box& yflxbx = amrex::surroundingNodes(bx, cdir);
    amrex::ParallelFor(yflxbx,
                       [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                         cns_riemann_y(i, j, k, nfabfy, slope, prims, cls);
                       });
    //

    // z-direction
    cdir = 2;
    const Box& zslpbx = amrex::grow(bx, cdir, 1);
    amrex::ParallelFor(zslpbx,
                       [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                         cns_slope_z(i, j, k, slope, prims, cls);
                       });
    const Box& zflxbx = amrex::surroundingNodes(bx, cdir);
    amrex::ParallelFor(zflxbx,
                       [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                         cns_riemann_z(i, j, k, nfabfz, slope, prims, cls);
                       });

    // don't have to do this, but we could
    // qeli.clear(); // don't need them anymore
    slopeeli.clear();


  }


  AMREX_GPU_DEVICE AMREX_FORCE_INLINE Real limiter(Real dlft,
                                                  Real drgt) noexcept {
    Real dcen = Real(0.5) * (dlft + drgt);
    Real dsgn = Math::copysign(Real(1.0), dcen);
    Real slop = Real(2.0) * min(Math::abs(dlft), Math::abs(drgt));
    Real dlim = (dlft * drgt >= Real(0.0)) ? slop : Real(0.0);
    return dsgn * min(dlim, Math::abs(dcen));
  }

  AMREX_GPU_DEVICE AMREX_FORCE_INLINE void cns_slope_x(
      int i, int j, int k, Array4<Real> const& dq, Array4<Real const> const& q,
      PROB::ProbClosures const& closures) noexcept {
    Real cspeed = sqrt(closures.gamma * closures.Rspec * q(i, j, k, QT)) + 2.e-40;
    Real dlft =
        Real(0.5) * (q(i, j, k, QPRES) - q(i - 1, j, k, QPRES)) / cspeed -
        Real(0.5) * q(i, j, k, QRHO) * (q(i, j, k, QU) - q(i - 1, j, k, QU));
    Real drgt =
        Real(0.5) * (q(i + 1, j, k, QPRES) - q(i, j, k, QPRES)) / cspeed -
        Real(0.5) * q(i, j, k, QRHO) * (q(i + 1, j, k, QU) - q(i, j, k, QU));
    Real d0 = limiter(dlft, drgt);

    Real cs2 = cspeed * cspeed;
    dlft = (q(i, j, k, QRHO) - q(i - 1, j, k, QRHO)) -
          (q(i, j, k, QPRES) - q(i - 1, j, k, QPRES)) / cs2;
    drgt = (q(i + 1, j, k, QRHO) - q(i, j, k, QRHO)) -
          (q(i + 1, j, k, QPRES) - q(i, j, k, QPRES)) / cs2;
    Real d1 = limiter(dlft, drgt);

    dlft = Real(0.5) * (q(i, j, k, QPRES) - q(i - 1, j, k, QPRES)) / cspeed +
          Real(0.5) * q(i, j, k, QRHO) * (q(i, j, k, QU) - q(i - 1, j, k, QU));
    drgt = Real(0.5) * (q(i + 1, j, k, QPRES) - q(i, j, k, QPRES)) / cspeed +
          Real(0.5) * q(i, j, k, QRHO) * (q(i + 1, j, k, QU) - q(i, j, k, QU));
    Real d2 = limiter(dlft, drgt);

    dlft = q(i, j, k, QV) - q(i - 1, j, k, QV);
    drgt = q(i + 1, j, k, QV) - q(i, j, k, QV);
    Real d3 = limiter(dlft, drgt);

    dlft = q(i, j, k, QW) - q(i - 1, j, k, QW);
    drgt = q(i + 1, j, k, QW) - q(i, j, k, QW);
    Real d4 = limiter(dlft, drgt);

    dq(i, j, k, 0) = d0;
    dq(i, j, k, 1) = d1;
    dq(i, j, k, 2) = d2;
    dq(i, j, k, 3) = d3;
    dq(i, j, k, 4) = d4;
  }

  AMREX_GPU_DEVICE AMREX_FORCE_INLINE void cns_slope_y(
      int i, int j, int k, Array4<Real> const& dq, Array4<Real const> const& q,
      PROB::ProbClosures const& closures) noexcept {
    Real cspeed = sqrt(closures.gamma * closures.Rspec * q(i, j, k, QT)) + 1.e-40;
    Real dlft =
        Real(0.5) * (q(i, j, k, QPRES) - q(i, j - 1, k, QPRES)) / cspeed -
        Real(0.5) * q(i, j, k, QRHO) * (q(i, j, k, QV) - q(i, j - 1, k, QV));
    Real drgt =
        Real(0.5) * (q(i, j + 1, k, QPRES) - q(i, j, k, QPRES)) / cspeed -
        Real(0.5) * q(i, j, k, QRHO) * (q(i, j + 1, k, QV) - q(i, j, k, QV));
    Real d0 = limiter(dlft, drgt);

    Real cs2 = cspeed * cspeed;
    dlft = (q(i, j, k, QRHO) - q(i, j - 1, k, QRHO)) -
          (q(i, j, k, QPRES) - q(i, j - 1, k, QPRES)) / cs2;
    drgt = (q(i, j + 1, k, QRHO) - q(i, j, k, QRHO)) -
          (q(i, j + 1, k, QPRES) - q(i, j, k, QPRES)) / cs2;
    Real d1 = limiter(dlft, drgt);

    dlft = Real(0.5) * (q(i, j, k, QPRES) - q(i, j - 1, k, QPRES)) / cspeed +
          Real(0.5) * q(i, j, k, QRHO) * (q(i, j, k, QV) - q(i, j - 1, k, QV));
    drgt = Real(0.5) * (q(i, j + 1, k, QPRES) - q(i, j, k, QPRES)) / cspeed +
          Real(0.5) * q(i, j, k, QRHO) * (q(i, j + 1, k, QV) - q(i, j, k, QV));
    Real d2 = limiter(dlft, drgt);

    dlft = q(i, j, k, QU) - q(i, j - 1, k, QU);
    drgt = q(i, j + 1, k, QU) - q(i, j, k, QU);
    Real d3 = limiter(dlft, drgt);

    dlft = q(i, j, k, QW) - q(i, j - 1, k, QW);
    drgt = q(i, j + 1, k, QW) - q(i, j, k, QW);
    Real d4 = limiter(dlft, drgt);

    dq(i, j, k, 0) = d0;
    dq(i, j, k, 1) = d1;
    dq(i, j, k, 2) = d2;
    dq(i, j, k, 3) = d3;
    dq(i, j, k, 4) = d4;
  }

  AMREX_GPU_DEVICE AMREX_FORCE_INLINE void cns_slope_z(
      int i, int j, int k, Array4<Real> const& dq, Array4<Real const> const& q,
      PROB::ProbClosures const& closures) noexcept {
    Real cspeed = sqrt(closures.gamma * closures.Rspec * q(i, j, k, QT)) + 1.e-40;
    Real dlft =
        Real(0.5) * (q(i, j, k, QPRES) - q(i, j, k - 1, QPRES)) / cspeed -
        Real(0.5) * q(i, j, k, QRHO) * (q(i, j, k, QW) - q(i, j, k - 1, QW));
    Real drgt =
        Real(0.5) * (q(i, j, k + 1, QPRES) - q(i, j, k, QPRES)) / cspeed -
        Real(0.5) * q(i, j, k, QRHO) * (q(i, j, k + 1, QW) - q(i, j, k, QW));
    Real d0 = limiter(dlft, drgt);

    Real cs2 = cspeed * cspeed;
    dlft = (q(i, j, k, QRHO) - q(i, j, k - 1, QRHO)) -
          (q(i, j, k, QPRES) - q(i, j, k - 1, QPRES)) / cs2;
    drgt = (q(i, j, k + 1, QRHO) - q(i, j, k, QRHO)) -
          (q(i, j, k + 1, QPRES) - q(i, j, k, QPRES)) / cs2;
    Real d1 = limiter(dlft, drgt);

    dlft = Real(0.5) * (q(i, j, k, QPRES) - q(i, j, k - 1, QPRES)) / cspeed +
          Real(0.5) * q(i, j, k, QRHO) * (q(i, j, k, QW) - q(i, j, k - 1, QW));
    drgt = Real(0.5) * (q(i, j, k + 1, QPRES) - q(i, j, k, QPRES)) / cspeed +
          Real(0.5) * q(i, j, k, QRHO) * (q(i, j, k + 1, QW) - q(i, j, k, QW));
    Real d2 = limiter(dlft, drgt);

    dlft = q(i, j, k, QU) - q(i, j, k - 1, QU);
    drgt = q(i, j, k + 1, QU) - q(i, j, k, QU);
    Real d3 = limiter(dlft, drgt);

    dlft = q(i, j, k, QV) - q(i, j, k - 1, QV);
    drgt = q(i, j, k + 1, QV) - q(i, j, k, QV);
    Real d4 = limiter(dlft, drgt);

    dq(i, j, k, 0) = d0;
    dq(i, j, k, 1) = d1;
    dq(i, j, k, 2) = d2;
    dq(i, j, k, 3) = d3;
    dq(i, j, k, 4) = d4;
  }

  AMREX_GPU_DEVICE AMREX_FORCE_INLINE void riemann_prob(
      const Real gamma, const Real smallp, const Real /*smallr*/, const Real rl,
      const Real ul, const Real pl, const Real ut1l, const Real ut2l,
      const Real rr, const Real ur, const Real pr, const Real ut1r,
      const Real ut2r, Real& flxrho, Real& flxu, Real& flxut, Real& flxutt,
      Real& flxe) noexcept {
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
    flxe = ugdnv * (Real(0.5) * rgdnv *
                        (ugdnv * ugdnv + utrans1 * utrans1 + utrans2 * utrans2) +
                    pgdnv / (gamma - Real(1.)) + pgdnv);
  }

  AMREX_GPU_DEVICE AMREX_FORCE_INLINE void cns_riemann_x(
      int i, int j, int k, Array4<Real> const& fx, Array4<Real const> const& dq,
      Array4<Real const> const& q, PROB::ProbClosures const& closures) noexcept {
    Real cspeed =
        sqrt(closures.gamma * closures.Rspec * q(i - 1, j, k, QT)) + 1.e-40;
    Real rl = q(i - 1, j, k, QRHO) +
              Real(0.5) * ((dq(i - 1, j, k, 0) + dq(i - 1, j, k, 2)) / cspeed +
                          dq(i - 1, j, k, 1));
    rl = max(rl, 1.0e-40);
    Real ul = q(i - 1, j, k, QU) +
              Real(0.5) * ((dq(i - 1, j, k, 2) - dq(i - 1, j, k, 0)) /
                          q(i - 1, j, k, QRHO));
    Real pl = q(i - 1, j, k, QPRES) +
              Real(0.5) * (dq(i - 1, j, k, 0) + dq(i - 1, j, k, 2)) * cspeed;
    pl = max(pl, 1.0e-40);
    Real ut1l = q(i - 1, j, k, QV) + Real(0.5) * dq(i - 1, j, k, 3);
    Real ut2l = q(i - 1, j, k, QW) + Real(0.5) * dq(i - 1, j, k, 4);

    cspeed = sqrt(closures.gamma * closures.Rspec * q(i, j, k, QT)) + 1.e-40;
    Real rr =
        q(i, j, k, QRHO) -
        Real(0.5) * ((dq(i, j, k, 0) + dq(i, j, k, 2)) / cspeed + dq(i, j, k, 1));
    rr = max(rr, 1.0e-40);
    Real ur = q(i, j, k, QU) -
              Real(0.5) * ((dq(i, j, k, 2) - dq(i, j, k, 0)) / q(i, j, k, QRHO));
    Real pr = q(i, j, k, QPRES) -
              Real(0.5) * (dq(i, j, k, 0) + dq(i, j, k, 2)) * cspeed;
    pr = max(pr, 1.0e-40);
    Real ut1r = q(i, j, k, QV) - Real(0.5) * dq(i, j, k, 3);
    Real ut2r = q(i, j, k, QW) - Real(0.5) * dq(i, j, k, 4);

    riemann_prob(closures.gamma, 1.0e-40, 1.0e-40, rl, ul, pl, ut1l, ut2l, rr, ur,
                pr, ut1r, ut2r, fx(i, j, k, URHO), fx(i, j, k, UMX),
                fx(i, j, k, UMY), fx(i, j, k, UMZ), fx(i, j, k, UET));
  }

  AMREX_GPU_DEVICE AMREX_FORCE_INLINE void cns_riemann_y(
      int i, int j, int k, Array4<Real> const& fy, Array4<Real const> const& dq,
      Array4<Real const> const& q, PROB::ProbClosures const& closures) noexcept {
    Real cspeed =
        sqrt(closures.gamma * closures.Rspec * q(i, j - 1, k, QT)) + 1.e-40;
    Real rl = q(i, j - 1, k, QRHO) +
              Real(0.5) * ((dq(i, j - 1, k, 0) + dq(i, j - 1, k, 2)) / cspeed +
                          dq(i, j - 1, k, 1));
    rl = max(rl, 1.0e-40);
    Real ul = q(i, j - 1, k, QV) +
              Real(0.5) * ((dq(i, j - 1, k, 2) - dq(i, j - 1, k, 0)) /
                          q(i, j - 1, k, QRHO));
    Real pl = q(i, j - 1, k, QPRES) +
              Real(0.5) * (dq(i, j - 1, k, 0) + dq(i, j - 1, k, 2)) * cspeed;
    pl = max(pl, 1.0e-40);
    Real ut1l = q(i, j - 1, k, QU) + Real(0.5) * dq(i, j - 1, k, 3);
    Real ut2l = q(i, j - 1, k, QW) + Real(0.5) * dq(i, j - 1, k, 4);

    cspeed = sqrt(closures.gamma * closures.Rspec * q(i, j, k, QT)) + 1.e-40;
    Real rr =
        q(i, j, k, QRHO) -
        Real(0.5) * ((dq(i, j, k, 0) + dq(i, j, k, 2)) / cspeed + dq(i, j, k, 1));
    rr = max(rr, 1.0e-40);
    Real ur = q(i, j, k, QV) -
              Real(0.5) * ((dq(i, j, k, 2) - dq(i, j, k, 0)) / q(i, j, k, QRHO));
    Real pr = q(i, j, k, QPRES) -
              Real(0.5) * (dq(i, j, k, 0) + dq(i, j, k, 2)) * cspeed;
    pr = max(pr, 1.0e-40);
    Real ut1r = q(i, j, k, QU) - Real(0.5) * dq(i, j, k, 3);
    Real ut2r = q(i, j, k, QW) - Real(0.5) * dq(i, j, k, 4);

    riemann_prob(closures.gamma, 1.0e-40, 1.0e-40, rl, ul, pl, ut1l, ut2l, rr, ur,
                pr, ut1r, ut2r, fy(i, j, k, URHO), fy(i, j, k, UMY),
                fy(i, j, k, UMX), fy(i, j, k, UMZ), fy(i, j, k, UET));
  }

  AMREX_GPU_DEVICE AMREX_FORCE_INLINE void cns_riemann_z(
      int i, int j, int k, Array4<Real> const& fz, Array4<Real const> const& dq,
      Array4<Real const> const& q, PROB::ProbClosures const& closures) noexcept {
    Real cspeed =
        sqrt(closures.gamma * closures.Rspec * q(i, j, k + 1, QT)) + 1.e-40;
    Real rl = q(i, j, k - 1, QRHO) +
              Real(0.5) * ((dq(i, j, k - 1, 0) + dq(i, j, k - 1, 2)) / cspeed +
                          dq(i, j, k - 1, 1));
    rl = max(rl, 1.0e-40);
    Real ul = q(i, j, k - 1, QW) +
              Real(0.5) * ((dq(i, j, k - 1, 2) - dq(i, j, k - 1, 0)) /
                          q(i, j, k - 1, QRHO));
    Real pl = q(i, j, k - 1, QPRES) +
              Real(0.5) * (dq(i, j, k - 1, 0) + dq(i, j, k - 1, 2)) * cspeed;
    pl = max(pl, 1.0e-40);
    Real ut1l = q(i, j, k - 1, QU) + Real(0.5) * dq(i, j, k - 1, 3);
    Real ut2l = q(i, j, k - 1, QV) + Real(0.5) * dq(i, j, k - 1, 4);

    cspeed = sqrt(closures.gamma * closures.Rspec * q(i, j, k, QT)) + 1.e-40;
    Real rr =
        q(i, j, k, QRHO) -
        Real(0.5) * ((dq(i, j, k, 0) + dq(i, j, k, 2)) / cspeed + dq(i, j, k, 1));
    rr = max(rr, 1.0e-40);
    Real ur = q(i, j, k, QW) -
              Real(0.5) * ((dq(i, j, k, 2) - dq(i, j, k, 0)) / q(i, j, k, QRHO));
    Real pr = q(i, j, k, QPRES) -
              Real(0.5) * (dq(i, j, k, 0) + dq(i, j, k, 2)) * cspeed;
    pr = max(pr, 1.0e-40);
    Real ut1r = q(i, j, k, QU) - Real(0.5) * dq(i, j, k, 3);
    Real ut2r = q(i, j, k, QV) - Real(0.5) * dq(i, j, k, 4);

    riemann_prob(closures.gamma, 1.0e-40, 1.0e-40, rl, ul, pl, ut1l, ut2l, rr, ur,
                pr, ut1r, ut2r, fz(i, j, k, URHO), fz(i, j, k, UMZ),
                fz(i, j, k, UMX), fz(i, j, k, UMY), fz(i, j, k, UET));
  }
}

#endif