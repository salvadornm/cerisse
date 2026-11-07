#ifndef ReconsHLLC_H_
#define ReconsHLLC_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_MFIter.H>

#include <cmath>

#include "recon.H"

/// \brief Reconstructed HLLC inviscid flux.
///
/// This class follows the Cerisse explicit RHS interface used by riemann_t/skew_t:
///
///   using ProbRHS = rhs_dt<reconshllc_t<param_t, ProbClosures>,
///                         no_diffusive_t, no_source_t>;
///
/// Required param members:
///   static constexpr int  recon_scheme = 1..6;  // see recon.H
///   static constexpr int  recon_sys    = 0 or 1; // 0: speed of sound, 1: gamma
///   static constexpr bool recon_char_var = true; // reconstruct characteristic vars
///   static constexpr Real plm_theta   = 1.5_rt; // used by MUSCL/PLM
///   static constexpr Real teno_cutoff = 0.0_rt; // used by TENO
///
/// Notes:
/// - If recon_char_var=false, primitive variables are reconstructed directly.
/// - With IBM/EB markers, fluxes at solid/covered faces are skipped, matching Riemann.h.
template <typename param, typename cls_t>
class reconshllc_t {
 public:
  AMREX_GPU_HOST_DEVICE
  reconshllc_t() {}

  AMREX_GPU_HOST_DEVICE
  ~reconshllc_t() {}

  static constexpr int recon_scheme = param::recon_scheme;
  static constexpr int recon_sys = param::recon_sys;
  static constexpr bool recon_char_var = param::recon_char_var;

  // Local characteristic-variable layout.  It mirrors hydro.H but avoids relying
  // on Pele index_macros.H names inside Cerisse problem files.
  static constexpr int WRHO_L  = 0;
  static constexpr int WACO_L  = 1; // WACO_L and WACO_L+1 are acoustic variables
  static constexpr int WTHER_L = 3;
  static constexpr int WUT_L   = 4;
  static constexpr int WY_L    = 4 + (AMREX_SPACEDIM - 1);
  static constexpr int NW      = WY_L + NUM_SPECIES;

#if (AMREX_USE_GPIBM || CNS_USE_EB)
  void inline eflux_ibm(const Geometry& /*geom*/, const MFIter& mfi,
                        const Array4<Real>& prims,
                        std::array<FArrayBox*, AMREX_SPACEDIM> const& flxt,
                        const Array4<Real>& /*rhs*/, const cls_t* cls,
                        const Array4<uint8_t>& ibMarkers)
#else
  void inline eflux(const Geometry& /*geom*/, const MFIter& mfi,
                    const Array4<Real>& prims,
                    std::array<FArrayBox*, AMREX_SPACEDIM> const& flxt,
                    const Array4<Real>& /*rhs*/, const cls_t* cls)
#endif
  {
    const Box& bx = mfi.tilebox();
    const int ngh = cls_t::NGHOST;

    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
      const Box flxbx = amrex::surroundingNodes(bx, dir);
      const Box recon_cell_bx = amrex::grow(bx, dir, 1);
      const Box face_store_bx = amrex::grow(flxbx, dir, 1);
      const Box char_bx = amrex::grow(bx, ngh);

      FArrayBox wf(char_bx, NW, The_Async_Arena());
      FArrayBox wlf(face_store_bx, NW, The_Async_Arena());
      FArrayBox wrf(face_store_bx, NW, The_Async_Arena());

      auto const& w = wf.array();
      auto const& wl = wlf.array();
      auto const& wr = wrf.array();
      auto const& flx = flxt[dir]->array();

      // Build either characteristic variables or primitive variables in w.
      amrex::ParallelFor(char_bx, [=, *this] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        if constexpr (recon_char_var) {
          this->cns_ctochar(i, j, k, dir, prims, w, *cls);
        } else {
          this->cns_ctoprim_recon(i, j, k, dir, prims, w, *cls);
        }
      });

      // Reconstruct left/right states at faces.
      amrex::ParallelFor(recon_cell_bx, NW,
        [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
          cns_recon<recon_scheme>(i, j, k, n, dir, w, wl, wr,
                                  param::plm_theta, param::teno_cutoff);
        });

      // HLLC flux on faces.
      amrex::ParallelFor(flxbx, [=, *this] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
#if (AMREX_USE_GPIBM || CNS_USE_EB)
        bool wallflx = false;
        if (dir == 0) {
          wallflx = ibMarkers(i, j, k, 0) && ibMarkers(i - 1, j, k, 0);
        }
#if AMREX_SPACEDIM >= 2
        else if (dir == 1) {
          wallflx = ibMarkers(i, j, k, 0) && ibMarkers(i, j - 1, k, 0);
        }
#endif
#if AMREX_SPACEDIM == 3
        else {
          wallflx = ibMarkers(i, j, k, 0) && ibMarkers(i, j, k - 1, 0);
        }
#endif
        if (!wallflx)
#endif
        {
          this->cns_riemann(i, j, k, dir, flx, prims, wl, wr, *cls);
        }
      });
    }
  }

  AMREX_GPU_DEVICE AMREX_FORCE_INLINE void cns_ctoprim_recon(
      int i, int j, int k, int dir, Array4<Real const> const& q,
      Array4<Real> const& w, cls_t const& cls) const noexcept {
    const int QU1 = cls.QU + dir;
    const int QU2 = (dir == 0) ? cls.QV : cls.QU;
    const int QU3 = (dir == 2) ? cls.QV : cls.QW;

    w(i, j, k, WRHO_L) = q(i, j, k, cls.QRHO);
    w(i, j, k, WACO_L) = q(i, j, k, QU1);
    w(i, j, k, WACO_L + 1) = q(i, j, k, cls.QPRES);
    w(i, j, k, WTHER_L) = q(i, j, k, cls.QC);
    AMREX_D_TERM(, w(i, j, k, WUT_L) = q(i, j, k, QU2);,
                   w(i, j, k, WUT_L + 1) = q(i, j, k, QU3););
    for (int n = 0; n < NUM_SPECIES; ++n) {
      w(i, j, k, WY_L + n) = q(i, j, k, cls.QFS + n);
    }
  }

  AMREX_GPU_DEVICE AMREX_FORCE_INLINE void cns_ctochar(
      int i, int j, int k, int dir, Array4<Real const> const& q,
      Array4<Real> const& w, cls_t const& cls) const noexcept {
    const int QU1 = cls.QU + dir;
    const int QU2 = (dir == 0) ? cls.QV : cls.QU;
    const int QU3 = (dir == 2) ? cls.QV : cls.QW;

    if constexpr (recon_sys == 0) {
      w(i, j, k, WRHO_L) = q(i, j, k, cls.QRHO) -
        q(i, j, k, cls.QPRES) / q(i, j, k, cls.QC) / q(i, j, k, cls.QC);
      w(i, j, k, WACO_L) = 0.5_rt * (q(i, j, k, cls.QPRES) / q(i, j, k, cls.QC) +
                                     q(i, j, k, cls.QRHO) * q(i, j, k, QU1));
      w(i, j, k, WACO_L + 1) = 0.5_rt * (q(i, j, k, cls.QPRES) / q(i, j, k, cls.QC) -
                                         q(i, j, k, cls.QRHO) * q(i, j, k, QU1));
      w(i, j, k, WTHER_L) = q(i, j, k, cls.QC);
    } else {
      w(i, j, k, WRHO_L) = q(i, j, k, cls.QRHO) * (1.0_rt - 1.0_rt / q(i, j, k, cls.QG));
      w(i, j, k, WACO_L) = 0.5_rt *
        (q(i, j, k, cls.QPRES) +
         std::sqrt(q(i, j, k, cls.QG) * q(i, j, k, cls.QRHO) * q(i, j, k, cls.QPRES)) *
           q(i, j, k, QU1));
      w(i, j, k, WACO_L + 1) = 0.5_rt *
        (q(i, j, k, cls.QPRES) -
         std::sqrt(q(i, j, k, cls.QG) * q(i, j, k, cls.QRHO) * q(i, j, k, cls.QPRES)) *
           q(i, j, k, QU1));
      w(i, j, k, WTHER_L) = q(i, j, k, cls.QG);
    }

    for (int n = 0; n < NUM_SPECIES; ++n) {
      w(i, j, k, WY_L + n) = q(i, j, k, cls.QFS + n);
    }
    AMREX_D_TERM(, w(i, j, k, WUT_L) = q(i, j, k, QU2);,
                   w(i, j, k, WUT_L + 1) = q(i, j, k, QU3););
  }

  AMREX_GPU_DEVICE AMREX_FORCE_INLINE void cns_unpackchar(
      Real w1, Real w2, Real w3, Real c, Real gamma,
      Real& rho, Real& u, Real& p) const noexcept {
    if constexpr (recon_sys == 0) {
      rho = w1 + (w2 + w3) / c;
      u = (w2 - w3) / rho;
      p = (w2 + w3) * c;
    } else {
      p = w2 + w3;
      rho = w1 / (1.0_rt - 1.0_rt / gamma);
      u = (w2 - w3) / std::sqrt(gamma * rho * p);
    }
  }

  AMREX_GPU_DEVICE AMREX_FORCE_INLINE void normalize_Y(Real Y[NUM_SPECIES]) const noexcept {
    Real sumY = 0.0_rt;
    for (int n = 0; n < NUM_SPECIES; ++n) {
      Y[n] = amrex::max(0.0_rt, Y[n]);
      sumY += Y[n];
    }
    if (sumY > 0.0_rt) {
      const Real invsumY = 1.0_rt / sumY;
      for (int n = 0; n < NUM_SPECIES; ++n) { Y[n] *= invsumY; }
    } else {
      Y[0] = 1.0_rt;
      for (int n = 1; n < NUM_SPECIES; ++n) { Y[n] = 0.0_rt; }
    }
  }

  AMREX_GPU_DEVICE AMREX_FORCE_INLINE void hllc(
      const Real rl, const Real ul, const Real pl, const Real ut1l,
      const Real ut2l, const Real el, const Real yl[NUM_SPECIES], const Real cl,
      const Real rr, const Real ur, const Real pr, const Real ut1r,
      const Real ut2r, const Real er, const Real yr[NUM_SPECIES], const Real cr,
      Real& flxrho, Real& flxu, Real& flxut, Real& flxutt, Real& flxrhoe,
      Real flxrhoy[NUM_SPECIES]) const noexcept {
    Real sl = amrex::min(ul - cl, ur - cr);
    Real sr = amrex::max(ul + cl, ur + cr);
    const Real rp = std::sqrt(rr / rl);
    const Real uroe = (ul + ur * rp) / (1.0_rt + rp);
    const Real croe = (cl + cr * rp) / (1.0_rt + rp);
    sl = amrex::min(sl, uroe - croe);
    sr = amrex::max(sr, uroe + croe);

    flxrho = 0.0_rt;
    if (sl > 0.0_rt) {
      flxu = rl * ul * ul + pl;
      flxut = rl * ul * ut1l;
      flxutt = rl * ul * ut2l;
      flxrhoe = ul * (rl * el + pl);
      for (int n = 0; n < NUM_SPECIES; ++n) {
        flxrhoy[n] = rl * ul * yl[n];
        flxrho += flxrhoy[n];
      }
    } else if (sr < 0.0_rt) {
      flxu = rr * ur * ur + pr;
      flxut = rr * ur * ut1r;
      flxutt = rr * ur * ut2r;
      flxrhoe = ur * (rr * er + pr);
      for (int n = 0; n < NUM_SPECIES; ++n) {
        flxrhoy[n] = rr * ur * yr[n];
        flxrho += flxrhoy[n];
      }
    } else {
      const Real sstar = (pr - pl + rl * ul * (sl - ul) - rr * ur * (sr - ur)) /
                         (rl * (sl - ul) - rr * (sr - ur));
      if (sstar >= 0.0_rt) {
        const Real frac = (sl - ul) / (sl - sstar) - 1.0_rt;
        flxu = rl * ul * ul + pl + sl * rl * ((frac + 1.0_rt) * sstar - ul);
        flxut = rl * ul * ut1l + sl * rl * frac * ut1l;
        flxutt = rl * ul * ut2l + sl * rl * frac * ut2l;
        flxrhoe = ul * (rl * el + pl) +
          sl * rl * (frac * el + (sstar - ul) * (sstar + pl / rl / (sl - ul)));
        for (int n = 0; n < NUM_SPECIES; ++n) {
          flxrhoy[n] = rl * ul * yl[n] + sl * rl * frac * yl[n];
          flxrho += flxrhoy[n];
        }
      } else {
        const Real frac = (sr - ur) / (sr - sstar) - 1.0_rt;
        flxu = rr * ur * ur + pr + sr * rr * ((frac + 1.0_rt) * sstar - ur);
        flxut = rr * ur * ut1r + sr * rr * frac * ut1r;
        flxutt = rr * ur * ut2r + sr * rr * frac * ut2r;
        flxrhoe = ur * (rr * er + pr) +
          sr * rr * (frac * er + (sstar - ur) * (sstar + pr / rr / (sr - ur)));
        for (int n = 0; n < NUM_SPECIES; ++n) {
          flxrhoy[n] = rr * ur * yr[n] + sr * rr * frac * yr[n];
          flxrho += flxrhoy[n];
        }
      }
    }
  }

  AMREX_GPU_DEVICE AMREX_FORCE_INLINE void cns_riemann(
      int i, int j, int k, int dir, Array4<Real> const& flx,
      Array4<Real const> const& q, Array4<Real> const& wl,
      Array4<Real> const& wr, cls_t const& cls) const noexcept {
    const amrex::IntVect iv(AMREX_D_DECL(i, j, k));
    const amrex::IntVect iv_dir = amrex::IntVect::TheDimensionVector(dir);

    const int QU1 = cls.QU + dir;
    const int QU2 = (dir == 0) ? cls.QV : cls.QU;
    const int QU3 = (dir == 2) ? cls.QV : cls.QW;
    const int UM1 = cls.UMX + dir;
    const int UM2 = (dir == 0) ? cls.UMY : cls.UMX;
    const int UM3 = (dir == 2) ? cls.UMY : cls.UMZ;

    Real rl, ul, pl, rr, ur, pr;
    Real ut1l = 0.0_rt, ut2l = 0.0_rt, ut1r = 0.0_rt, ut2r = 0.0_rt;
    Real yl[NUM_SPECIES], yr[NUM_SPECIES];

    if constexpr (recon_char_var) {
      if (!(wl(iv, WTHER_L) > 0.0_rt) || !(wr(iv, WTHER_L) > 0.0_rt)) {
        wl(iv, WTHER_L) = (recon_sys == 0) ? q(iv - iv_dir, cls.QC) : q(iv - iv_dir, cls.QG);
        wr(iv, WTHER_L) = (recon_sys == 0) ? q(iv, cls.QC) : q(iv, cls.QG);
      }
      const Real gl = (recon_sys == 1) ? wl(iv, WTHER_L) : q(iv - iv_dir, cls.QG);
      const Real cl0 = (recon_sys == 0) ? wl(iv, WTHER_L) : q(iv - iv_dir, cls.QC);
      cns_unpackchar(wl(iv, WRHO_L), wl(iv, WACO_L), wl(iv, WACO_L + 1),
                     cl0, gl, rl, ul, pl);

      const Real gr = (recon_sys == 1) ? wr(iv, WTHER_L) : q(iv, cls.QG);
      const Real cr0 = (recon_sys == 0) ? wr(iv, WTHER_L) : q(iv, cls.QC);
      cns_unpackchar(wr(iv, WRHO_L), wr(iv, WACO_L), wr(iv, WACO_L + 1),
                     cr0, gr, rr, ur, pr);

      AMREX_D_TERM(, ut1l = wl(iv, WUT_L); ut1r = wr(iv, WUT_L);,
                     ut2l = wl(iv, WUT_L + 1); ut2r = wr(iv, WUT_L + 1););
      for (int n = 0; n < NUM_SPECIES; ++n) {
        yl[n] = wl(iv, WY_L + n);
        yr[n] = wr(iv, WY_L + n);
      }
    } else {
      rl = wl(iv, WRHO_L);  ul = wl(iv, WACO_L);      pl = wl(iv, WACO_L + 1);
      rr = wr(iv, WRHO_L);  ur = wr(iv, WACO_L);      pr = wr(iv, WACO_L + 1);
      AMREX_D_TERM(, ut1l = wl(iv, WUT_L); ut1r = wr(iv, WUT_L);,
                     ut2l = wl(iv, WUT_L + 1); ut2r = wr(iv, WUT_L + 1););
      for (int n = 0; n < NUM_SPECIES; ++n) {
        yl[n] = wl(iv, WY_L + n);
        yr[n] = wr(iv, WY_L + n);
      }
    }

    normalize_Y(yl);
    normalize_Y(yr);

    if (!(rl > 0.0_rt) || !(pl > 0.0_rt) || !(rr > 0.0_rt) || !(pr > 0.0_rt) ||
        amrex::isnan(AMREX_D_TERM(ul + ur, +ut1l + ut1r, +ut2l + ut2r))) {
      rl = q(iv - iv_dir, cls.QRHO);
      ul = q(iv - iv_dir, QU1);
      pl = q(iv - iv_dir, cls.QPRES);
      rr = q(iv, cls.QRHO);
      ur = q(iv, QU1);
      pr = q(iv, cls.QPRES);
      

      AMREX_D_TERM(, ut1l = q(iv - iv_dir, QU2); ut1r = q(iv, QU2);,
                     ut2l = q(iv - iv_dir, QU3); ut2r = q(iv, QU3););
      for (int n = 0; n < NUM_SPECIES; ++n) {
        yl[n] = q(iv - iv_dir, cls.QFS + n);
        yr[n] = q(iv, cls.QFS + n);
      }
#if NUM_SPECIES > 1      
      normalize_Y(yl);
      normalize_Y(yr);
#endif      
    }

    Real el, er, cl, cr;
    cls.RYP2E(rl, yl, pl, el);
    el += 0.5_rt * (AMREX_D_TERM(ul * ul, +ut1l * ut1l, +ut2l * ut2l));
    cl = q(iv - iv_dir, cls.QC);
    

    cls.RYP2E(rr, yr, pr, er);
    er += 0.5_rt * (AMREX_D_TERM(ur * ur, +ut1r * ut1r, +ut2r * ut2r));
    cr = q(iv, cls.QC);

    // if (amrex::isnan(cl) || amrex::isnan(cr)) {
    //   cl = q(iv - iv_dir, cls.QC);
    //   cr = q(iv, cls.QC);
    // }

    Real flxrY[NUM_SPECIES] = {0.0_rt};
    hllc(rl, ul, pl, ut1l, ut2l, el, yl, cl,
         rr, ur, pr, ut1r, ut2r, er, yr, cr,
         flx(iv, cls.URHO), flx(iv, UM1), flx(iv, UM2), flx(iv, UM3),
         flx(iv, cls.UET), flxrY);

    for (int n = 0; n < NUM_SPECIES; ++n) {
      flx(iv, cls.UFS + n) = flxrY[n];
    }
  }

  // self-consistent HLLC piecewise
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE void cns_riemann_first_order(
    amrex::IntVect const& iv_face, int dir,
    amrex::Array4<amrex::Real> const& flx, amrex::Real const* ql, amrex::Real const* qr,
    cls_t const& cls) const noexcept
    {
      const int QU1 = cls.QU  + dir;
      const int UM1 = cls.UMX + dir;

      const int QU2 = (dir == 0) ? cls.QV : cls.QU;
      const int UM2 = (dir == 0) ? cls.UMY : cls.UMX;

      const int QU3 = (dir == 2) ? cls.QV : cls.QW;
      const int UM3 = (dir == 2) ? cls.UMY : cls.UMZ;

      const Real rl = ql[cls.QRHO]; const Real ul = ql[QU1]; const Real pl = ql[cls.QPRES]; const Real cl = ql[cls.QC];

      const Real rr = qr[cls.QRHO]; const Real ur = qr[QU1]; const Real pr = qr[cls.QPRES]; const Real cr = qr[cls.QC];

      Real ut1l = 0.0_rt; Real ut1r = 0.0_rt; Real ut2l = 0.0_rt;Real ut2r = 0.0_rt;

#if AMREX_SPACEDIM >= 2
      ut1l = ql[QU2];
      ut1r = qr[QU2];
#endif
#if AMREX_SPACEDIM == 3
      ut2l = ql[QU3];
      ut2r = qr[QU3];
#endif

      Real Yl[NUM_SPECIES];
      Real Yr[NUM_SPECIES];

      for (int ns = 0; ns < NUM_SPECIES; ++ns) {
        Yl[ns] = ql[cls.QFS + ns];
        Yr[ns] = qr[cls.QFS + ns];
      }

#if NUM_SPECIES > 1
      normalize_Y(Yl);
      normalize_Y(Yr);
#endif

      Real el, er;
      cls.RYP2E(rl, Yl, pl, el);
      cls.RYP2E(rr, Yr, pr, er);

      el += Real(0.5) *
          (AMREX_D_TERM(ul*ul, +ut1l*ut1l, +ut2l*ut2l));

      er += Real(0.5) *
          (AMREX_D_TERM(ur*ur, +ut1r*ut1r, +ut2r*ut2r));

      Real fluxY[NUM_SPECIES] = {Real(0.0)};

      hllc(
        rl, ul, pl, ut1l, ut2l, el, Yl, cl,
        rr, ur, pr, ut1r, ut2r, er, Yr, cr,
        flx(iv_face, cls.URHO),
        flx(iv_face, UM1),
        flx(iv_face, UM2),
        flx(iv_face, UM3),
        flx(iv_face, cls.UET),
        fluxY);

      for (int ns = 0; ns < NUM_SPECIES; ++ns) {
        flx(iv_face, cls.UFS + ns) = fluxY[ns];
      }
    }
  //
};

// Optional alias with the spelling used in the question.
template <typename param, typename cls_t>
using reconHLLC_t = reconshllc_t<param, cls_t>;

template <typename param, typename cls_t>
using reconshllc_t_t = reconshllc_t<param, cls_t>;

#endif
