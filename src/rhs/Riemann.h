#ifndef Riemann_H_
#define Riemann_H_

#include <AMReX_FArrayBox.H>
#include <CNS.h>

///
/// \brief Template class for Riemann solvers
template <bool iOption, typename cls_t>
class riemann_t {
 public:
  AMREX_GPU_HOST_DEVICE
  riemann_t() {}

  AMREX_GPU_HOST_DEVICE
  ~riemann_t() {}

#if (AMREX_USE_GPIBM || CNS_USE_EB )  
 void inline eflux_ibm(const Geometry& geom, const MFIter& mfi,
                    const Array4<Real>& prims, std::array<FArrayBox*, AMREX_SPACEDIM> const &flxt,
                    const Array4<Real>& rhs, const cls_t* cls,const Array4<bool>& ibMarkers) {

#else
  void inline eflux(const Geometry& geom, const MFIter& mfi,
                    const Array4<Real>& prims, std::array<FArrayBox*, AMREX_SPACEDIM> const &flxt,
                    const Array4<Real>& rhs, const cls_t* cls) {
#endif
  
    const Box& bx  = mfi.tilebox();
    const Box& bxg = mfi.growntilebox(cls_t::NGHOST);

    // zero rhs
    ParallelFor(bxg, cls_t::NCONS,
                [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
                  rhs(i, j, k, n) = 0.0;
                });
    const Box& bxg1 = amrex::grow(bx, 1);
    FArrayBox slopef(bxg1, cls_t::NCONS+1, The_Async_Arena()); // +1 because of the density
    const Array4<Real>& slope = slopef.array();

    // x-direction
    int cdir = 0;
    auto const& flx1 = flxt[cdir]->array(); 
  

    const Box& xslpbx = amrex::grow(bx, cdir, 1);


    amrex::ParallelFor(
        xslpbx, [=, *this] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
#ifdef CNS_USE_EB 
          this->cns_slope_x(i, j, k, slope, prims, *cls, ibMarkers); 
#else          
          this->cns_slope_x(i, j, k, slope, prims, *cls);
#endif          

        });


    const Box& xflxbx = amrex::surroundingNodes(bx, cdir);
    amrex::ParallelFor(
        xflxbx, [=, *this] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
          this->cns_riemann_x(i, j, k, flx1, slope, prims, *cls);
        });

#if AMREX_SPACEDIM >= 2
    // y-direction
    cdir = 1;
    auto const& flx2 = flxt[cdir]->array(); 
  
    const Box& yslpbx = amrex::grow(bx, cdir, 1);
    amrex::ParallelFor(
        yslpbx, [=, *this] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {        
#ifdef CNS_USE_EB 
          this->cns_slope_y(i, j, k, slope, prims, *cls, ibMarkers); 
#else          
          this->cns_slope_y(i, j, k, slope, prims, *cls);
#endif  

        });
    const Box& yflxbx = amrex::surroundingNodes(bx, cdir);
    amrex::ParallelFor(
        yflxbx, [=, *this] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {          
          this->cns_riemann_y(i, j, k, flx2, slope, prims, *cls);
        });
#endif

#if AMREX_SPACEDIM == 3
    // z-direction
    cdir = 2;
    auto const& flx3 = flxt[cdir]->array(); 
  
    const Box& zslpbx = amrex::grow(bx, cdir, 1);
    amrex::ParallelFor(
        zslpbx, [=, *this] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
#ifdef CNS_USE_EB 
          this->cns_slope_z(i, j, k, slope, prims, *cls, ibMarkers); 
#else          
          this->cns_slope_z(i, j, k, slope, prims, *cls);
#endif 


        });
    const Box& zflxbx = amrex::surroundingNodes(bx, cdir);
    amrex::ParallelFor(
        zflxbx, [=, *this] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
          this->cns_riemann_z(i, j, k, flx3, slope, prims, *cls);
        });
#endif
  };

  AMREX_GPU_DEVICE AMREX_FORCE_INLINE Real limiter(Real dlft, Real drgt) const {
    Real dcen = Real(0.5) * (dlft + drgt);
    Real dsgn = Math::copysign(Real(1.0), dcen);
    Real slop = Real(2.0) * min(Math::abs(dlft), Math::abs(drgt));
    Real dlim = (dlft * drgt >= Real(0.0)) ? slop : Real(0.0);
    return dsgn * min(dlim, Math::abs(dcen));
  }

  AMREX_GPU_DEVICE AMREX_FORCE_INLINE void hllc(
      const Real rl, const Real ul, const Real pl, const Real ut1l,
      const Real ut2l, const Real el, const Real yl[NUM_SPECIES], const Real cl, 
      const Real rr, const Real ur, const Real pr, const Real ut1r, 
      const Real ut2r, const Real er, const Real yr[NUM_SPECIES], const Real cr, 
      Real& flxu, Real& flxut, Real& flxutt, Real& flxrhoe, Real flxrhoy[NUM_SPECIES]) const {

    // HLLC borrowed from multispecies branch

    using amrex::Real;

    // Estimate wave speeds
    Real sl = amrex::min(ul - cl, ur - cr);
    Real sr = amrex::max(ul + cl, ur + cr);
    Real rp = std::sqrt(rr / rl);
    Real uroe = (ul + ur * rp) / (1. + rp);
    Real croe = (cl + cr * rp) / (1. + rp); // a more standard averaging
    sl = amrex::min(sl, uroe - croe);
    sr = amrex::max(sr, uroe + croe);

    if (sl > 0) {
      // flx_l
      flxu = rl * ul * ul + pl;
      flxut = rl * ul * ut1l;
      flxutt = rl * ul * ut2l;
      flxrhoe = ul * (rl * el + pl);
      for (int n = 0; n < NUM_SPECIES; ++n) {
        flxrhoy[n] = rl * ul * yl[n];
      }

    } else if (sr < 0) {
      // flx_r
      flxu = rr * ur * ur + pr;
      flxut = rr * ur * ut1r;
      flxutt = rr * ur * ut2r;
      flxrhoe = ur * (rr * er + pr);
      for (int n = 0; n < NUM_SPECIES; ++n) {
        flxrhoy[n] = rr * ur * yr[n];
      }

    } else {
      Real sstar = (pr - pl + rl * ul * (sl - ul) - rr * ur * (sr - ur)) /
                  (rl * (sl - ul) - rr * (sr - ur)); // contact wave speed

      if (sstar >= 0) {
        // flx_l* = flx_l + sl * (q_l* - q_l)
        Real frac = (sl - ul) / (sl - sstar) - 1.;

        flxu = rl * ul * ul + pl + sl * rl * ((frac + 1.) * sstar - ul);
        flxut = rl * ul * ut1l + sl * rl * frac * ut1l;
        flxutt = rl * ul * ut2l + sl * rl * frac * ut2l;
        flxrhoe = ul * (rl * el + pl) +
                  sl * rl * (frac * el + (sstar - ul) * (sstar + pl / rl / (sl - ul)));
        for (int n = 0; n < NUM_SPECIES; ++n) {
          flxrhoy[n] = rl * ul * yl[n] + sl * rl * frac * yl[n];
        }

      } else {
        // flx_r* = flx_r + sr * (q_r* - q_r)
        Real frac = (sr - ur) / (sr - sstar) - 1.;

        flxu = rr * ur * ur + pr + sr * rr * ((frac + 1.) * sstar - ur);
        flxut = rr * ur * ut1r + sr * rr * frac * ut1r;
        flxutt = rr * ur * ut2r + sr * rr * frac * ut2r;
        flxrhoe = ur * (rr * er + pr) +
                  sr * rr * (frac * er + (sstar - ur) * (sstar + pr / rr / (sr - ur)));
        for (int n = 0; n < NUM_SPECIES; ++n) {
          flxrhoy[n] = rr * ur * yr[n] + sr * rr * frac * yr[n];
        }
      }
    }
  }

#ifdef CNS_USE_EB 
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE void cns_slope_x(
      int i, int j, int k, const Array4<Real>& dq, const Array4<Real>& q,
      const cls_t& cls, const Array4<bool>& marker) const {
#else
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE void cns_slope_x(
      int i, int j, int k, const Array4<Real>& dq, const Array4<Real>& q,
      const cls_t& cls) const {
#endif




#ifdef CNS_USE_EB 
    // const bool wall_left  = marker(i-1,j,k,0);  
    // const bool wall_right = marker(i+1,j,k,0);  
    const bool close_to_wall = marker(i,j,k,1);
#endif

    Real cspeed = q(i, j, k, cls.QC) + 1.e-40;
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

    for (int n = 0; n < NUM_SPECIES; ++n) {
      dlft = q(i, j, k, cls.QFS + n) - q(i - 1, j, k, cls.QFS + n);
      drgt = q(i + 1, j, k, cls.QFS + n) - q(i, j, k, cls.QFS + n);
      dq(i, j, k, 5 + n) = limiter(dlft, drgt);
    }

#ifdef CNS_USE_EB
    if (close_to_wall) 
    {
      for (int n = 0; n < cls.NCONS; ++n) { dq(i,j,k,n) = 0.0;}
    }    
#endif    

  }
 
#ifdef CNS_USE_EB 
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE void cns_slope_y(
      int i, int j, int k, const Array4<Real>& dq, const Array4<Real>& q,
      const cls_t& cls, const Array4<bool>& marker) const {
#else
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE void cns_slope_y(
      int i, int j, int k, const Array4<Real>& dq, const Array4<Real>& q,
      const cls_t& cls) const {
#endif



#ifdef CNS_USE_EB 
    // const bool wall_left  = marker(i,j-1,k,0);  
    // const bool wall_right = marker(i,j+1,k,0); 
    const bool close_to_wall = marker(i,j,k,1); 
#endif


    Real cspeed = q(i, j, k, cls.QC) + 1.e-40;
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

    for (int n = 0; n < NUM_SPECIES; ++n) {
      dlft = q(i, j, k, cls.QFS + n) - q(i, j - 1, k, cls.QFS + n);
      drgt = q(i, j + 1, k, cls.QFS + n) - q(i, j, k, cls.QFS + n);
      dq(i, j, k, 5 + n) = limiter(dlft, drgt);
    }

#ifdef CNS_USE_EB
    if (close_to_wall) 
    {
      for (int n = 0; n < cls.NCONS; ++n) { dq(i,j,k,n) = 0.0;}
    }    
#endif    
 

  }

#ifdef CNS_USE_EB 
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE void cns_slope_z(
      int i, int j, int k, const Array4<Real>& dq, const Array4<Real>& q,
      const cls_t& cls, const Array4<bool>& marker) const {
#else
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE void cns_slope_z(
      int i, int j, int k, const Array4<Real>& dq, const Array4<Real>& q,
      const cls_t& cls) const {
#endif
  
#ifdef CNS_USE_EB 
    // const bool wall_left  = marker(i,j-1,k,0);  
    // const bool wall_right = marker(i,j+1,k,0); 
    const bool close_to_wall = marker(i,j,k,1); 
#endif

    Real cspeed = q(i, j, k, cls.QC) + 1.e-40;
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

    for (int n = 0; n < NUM_SPECIES; ++n) {
      dlft = q(i, j, k, cls.QFS + n) - q(i, j, k - 1, cls.QFS + n);
      drgt = q(i, j, k + 1, cls.QFS + n) - q(i, j, k, cls.QFS + n);
      dq(i, j, k, 5 + n) = limiter(dlft, drgt);
    }

#ifdef CNS_USE_EB
    if (close_to_wall) 
    {
      for (int n = 0; n < cls.NCONS; ++n) { dq(i,j,k,n) = 0.0;}
    }    
#endif   

  }

  AMREX_GPU_DEVICE AMREX_FORCE_INLINE void cns_riemann_x(
      int i, int j, int k, Array4<Real> const& fx, Array4<Real> const& dq,
      Array4<Real> const& q, cls_t const& cls) const {
    Real cspeed = q(i - 1, j, k, cls.QC) + 1.e-40;
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

    cspeed = q(i, j, k, cls.QC) + 1.e-40;
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

    Real Yl[NUM_SPECIES], Yr[NUM_SPECIES], flxrY[NUM_SPECIES];
    for (int n = 0; n < NUM_SPECIES; ++n) {
      Yl[n] = q(i - 1, j, k, cls.QFS + n) + Real(0.5) * dq(i - 1, j, k, 5 + n);
      Yr[n] = q(i, j, k, cls.QFS + n) - Real(0.5) * dq(i, j, k, 5 + n);
    }

    // const Real gamma = 0.5 * (q(i - 1, j, k, cls.QG) + q(i, j, k, cls.QG));
    // riemann_prob(gamma, 1.0e-40, 1.0e-40, rl, ul, pl, ut1l, ut2l, Yl, rr, ur,
    //              pr, ut1r, ut2r, Yr, fx(i, j, k, cls.UMX), fx(i, j, k, cls.UMY),
    //              fx(i, j, k, cls.UMZ), fx(i, j, k, cls.UET), flxrY);
    
    Real el, er;
    cls.RYP2E(rl, Yl, pl, el);
    el += 0.5 * (AMREX_D_TERM(ul * ul, +ut1l * ut1l, +ut2l * ut2l));
    cls.RYP2E(rr, Yr, pr, er);
    er += 0.5 * (AMREX_D_TERM(ur * ur, +ut1r * ut1r, +ut2r * ut2r));
    hllc(
      rl, ul, pl, ut1l, ut2l, el, Yl, q(i-1, j, k, cls.QC), 
      rr, ur, pr, ut1r, ut2r, er, Yr, q(i, j, k, cls.QC), 
      fx(i, j, k, cls.UMX), fx(i, j, k, cls.UMY),
      fx(i, j, k, cls.UMZ), fx(i, j, k, cls.UET), flxrY);

    for (int n = 0; n < NUM_SPECIES; ++n) {
      fx(i, j, k, cls.UFS + n) = flxrY[n];
    }
  }

  AMREX_GPU_DEVICE AMREX_FORCE_INLINE void cns_riemann_y(
      int i, int j, int k, Array4<Real> const& fy, Array4<Real const> const& dq,
      Array4<Real const> const& q, cls_t const& cls) const {
    Real cspeed = q(i, j - 1, k, cls.QC) + 1.e-40;
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

    cspeed = q(i, j, k, cls.QC) + 1.e-40;
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

    Real Yl[NUM_SPECIES], Yr[NUM_SPECIES], flxrY[NUM_SPECIES];
    for (int n = 0; n < NUM_SPECIES; ++n) {
      Yl[n] = q(i, j - 1, k, cls.QFS + n) + Real(0.5) * dq(i, j - 1, k, 5 + n);
      Yr[n] = q(i, j, k, cls.QFS + n) - Real(0.5) * dq(i, j, k, 5 + n);
    }
    // const Real gamma = 0.5 * (q(i, j - 1, k, cls.QG) + q(i, j, k, cls.QG));
    // riemann_prob(gamma, 1.0e-40, 1.0e-40, rl, ul, pl, ut1l, ut2l, Yl, rr, ur,
    //              pr, ut1r, ut2r, Yr, fy(i, j, k, cls.UMY), fy(i, j, k, cls.UMX),
    //              fy(i, j, k, cls.UMZ), fy(i, j, k, cls.UET), flxrY);
    
    Real el, er;
    cls.RYP2E(rl, Yl, pl, el);
    el += 0.5 * (AMREX_D_TERM(ul * ul, +ut1l * ut1l, +ut2l * ut2l));
    cls.RYP2E(rr, Yr, pr, er);
    er += 0.5 * (AMREX_D_TERM(ur * ur, +ut1r * ut1r, +ut2r * ut2r));
    hllc(
      rl, ul, pl, ut1l, ut2l, el, Yl, q(i, j-1, k, cls.QC), 
      rr, ur, pr, ut1r, ut2r, er, Yr, q(i, j, k, cls.QC), 
      fy(i, j, k, cls.UMY), fy(i, j, k, cls.UMX),
      fy(i, j, k, cls.UMZ), fy(i, j, k, cls.UET), flxrY);

    for (int n = 0; n < NUM_SPECIES; ++n) {
      fy(i, j, k, cls.UFS + n) = flxrY[n];
    }
  }

  AMREX_GPU_DEVICE AMREX_FORCE_INLINE void cns_riemann_z(
      int i, int j, int k, Array4<Real> const& fz, Array4<Real const> const& dq,
      Array4<Real const> const& q, cls_t const& cls) const {
    Real cspeed = q(i, j, k - 1, cls.QC) + 1.e-40;
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

    cspeed = q(i, j, k, cls.QC) + 1.e-40;
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

    Real Yl[NUM_SPECIES], Yr[NUM_SPECIES], flxrY[NUM_SPECIES];
    for (int n = 0; n < NUM_SPECIES; ++n) {
      Yl[n] = q(i, j, k - 1, cls.QFS + n) + Real(0.5) * dq(i, j, k - 1, 5 + n);
      Yr[n] = q(i, j, k, cls.QFS + n) - Real(0.5) * dq(i, j, k, 5 + n);
    }
    // const Real gamma = 0.5 * (q(i - 1, j, k, cls.QG) + q(i, j, k, cls.QG));
    // riemann_prob(gamma, 1.0e-40, 1.0e-40, rl, ul, pl, ut1l, ut2l, Yl, rr, ur,
    //              pr, ut1r, ut2r, Yr, fz(i, j, k, cls.UMZ), fz(i, j, k, cls.UMX),
    //              fz(i, j, k, cls.UMY), fz(i, j, k, cls.UET), flxrY);
    
    Real el, er;
    cls.RYP2E(rl, Yl, pl, el);
    el += 0.5 * (AMREX_D_TERM(ul * ul, +ut1l * ut1l, +ut2l * ut2l));
    cls.RYP2E(rr, Yr, pr, er);
    er += 0.5 * (AMREX_D_TERM(ur * ur, +ut1r * ut1r, +ut2r * ut2r));
    hllc(
      rl, ul, pl, ut1l, ut2l, el, Yl, q(i, j, k-1, cls.QC), 
      rr, ur, pr, ut1r, ut2r, er, Yr, q(i, j, k, cls.QC), 
      fz(i, j, k, cls.UMZ), fz(i, j, k, cls.UMX),
      fz(i, j, k, cls.UMY), fz(i, j, k, cls.UET), flxrY);

    for (int n = 0; n < NUM_SPECIES; ++n) {
      fz(i, j, k, cls.UFS + n) = flxrY[n];
    }
  }
  };

#endif