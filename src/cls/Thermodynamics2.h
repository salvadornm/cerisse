#ifndef THERMODYNAMICS2_H_
#define THERMODYNAMICS2_H_

using namespace amrex;

#include <Constants.h>
#include <CNSconstants.h>


using namespace universal_constants;
using namespace CNSConstants;

// Generic Perfect gas (needs gamma and moecular weight as input)
////////////////////////////////////////////////////////////////////////////////////////
template <typename param, typename idx_t>
class perfect_gas_t {
 protected:
 public:
  Real gamma   = param::gamma;        // ratio of specific heats
  Real mw = param::molecular_weight;  // mean molecular weight air kg/mol

  Real gamma_m1 = gamma - Real(1.0);
  Real o_gamma_m1 = Real(1.0)/gamma_m1;
  Real Ru = gas_constant;
  Real cv = Ru / (mw * gamma_m1);
  Real cp = gamma * cv;
  Real Rspec = Ru / mw;

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real get_ei_min() const {
#if CLIP_TEMPERATURE_MIN
    return Rspec * min_temp() * o_gamma_m1;
#else
    return min_press() * o_gamma_m1;
#endif
  }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void RYP2E(const Real R,
                                                      const Real* /*Y*/,
                                                      const Real P,
                                                      Real& E) const {
    E = P / (R * gamma_m1);
  }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void RYE2TP( const Real R,
                                                        const Real* /*Y*/,
                                                        const Real E, Real& T,
                                                        Real& P) const {
    P = gamma_m1 * R * E;
    T = P / (R * Rspec);
  }
  // \brief calculate the speed of sound (function of R,Y,E)
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void RYE2Cs( const Real /*R*/,
                                                        const Real* /*Y*/,
                                                        const Real E,
                                                        Real& cs) const {
    Real T = E / cv;
    cs = std::sqrt(gamma * Rspec * T);
  }
  // \brief calculate the speed of sound (function of P,YT)
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void PYT2Cs( const Real /*R*/,
                                                        const Real* /*Y*/,
                                                        const Real T,
                                                        Real& cs) const {
    cs = std::sqrt(gamma * Rspec * T);
  }


  // \brief calculate the density
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void PYT2R(
    const Real P, const Real* /*Y*/, const Real T, Real& R) const {

    R = P/(T *Rspec);    
  }    
  
  // \brief calculate the specific internal energy
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void PYT2E(
    const Real /*P*/, const Real* /*Y*/, const Real T, Real& E) const {      
    E = cv*T;
  }
  
  // \brief this function ensures P and T  do not violate bounds
  // and then fills the q-array  to ensure consistency.
  // Used in IBM to calculate auxiliar primitives corerctly

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void ensurePTYfillq(
    Real& P, Real& T, Real* /*Y*/, 
    const Real ux, const Real uy, const Real uz,
    Real* Q ) const {

    P = amrex::max(min_press(),P);  
#if CLIP_TEMPERATURE_MIN        
    T = amrex::max(min_temp(),T);
#endif    
    Q[idx_t::QRHO] = P/(T*Rspec);
    Q[idx_t::QT] = T;
    Q[idx_t::QPRES] = P;
    Q[idx_t::QU] = ux;
    Q[idx_t::QV] = uy;
    Q[idx_t::QW] = uz;
    Q[idx_t::QFS] = 1.0;
    // aux primitives
    Q[idx_t::QC] =  std::sqrt(gamma * Rspec * T);
    Q[idx_t::QG] = gamma; 
    Q[idx_t::QEINT] = cv*T;
    
  }   
  ////////////////////////////////////////////////////////////////

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE GpuArray<Real, idx_t::NWAVES>
  cons2eigenvals(const int i, const int j, const int k,
                 const Array4<Real>& cons, const GpuArray<int, 3>& vdir) const {

    Real rho = cons(i, j, k, idx_t::URHO);
    Real rhoinv = Real(1.0) / rho;
    GpuArray<Real, AMREX_SPACEDIM> vel = {AMREX_D_DECL(
        cons(i, j, k, idx_t::UMX) * rhoinv, cons(i, j, k, idx_t::UMY) * rhoinv,
        cons(i, j, k, idx_t::UMZ) * rhoinv)};

    Real ke = AMREX_D_PICK(vel[0] * vel[0],  vel[0] * vel[0] + vel[1] * vel[1],
                           vel[0] * vel[0] + vel[1] * vel[1] + vel[2] * vel[2]);
    ke = Real(0.5) * rho * ke;
    Real eint = (cons(i, j, k, idx_t::UET) - ke) / rho;
    eint = max(eint,get_ei_min() ); //clip energy
    Real T = eint / cv;

    Real cs = std::sqrt(gamma * Rspec * T);
    Real u =
        AMREX_D_PICK(vel[0] * vdir[0], vel[0] * vdir[0] + vel[1] * vdir[1],
                     vel[0] * vdir[0] + vel[1] * vdir[1] + vel[2] * vdir[2]);
    GpuArray<Real, idx_t::NWAVES> eigenvals = {u + cs, u, u - cs};
    return eigenvals;
  }

  AMREX_GPU_DEVICE AMREX_FORCE_INLINE void prims2fluxes(
      int i, int j, int k, const Array4<Real>& prims,
      const Array4<Real>& fluxes, const GpuArray<int, 3>& vdir) const {

    Real rho = prims(i, j, k, idx_t::QRHO);
    Real ux = prims(i, j, k, idx_t::QU);
    Real uy = prims(i, j, k, idx_t::QV);
    Real uz = prims(i, j, k, idx_t::QW);
    Real P = prims(i, j, k, idx_t::QPRES);
    Real udir = ux * vdir[0] + uy * vdir[1] + uz * vdir[2];

    Real ekin = Real(0.5) * (ux * ux + uy * uy + uz * uz);
    Real rhoet = rho * (cp * prims(i, j, k, idx_t::QT) + ekin);

    fluxes(i, j, k, idx_t::URHO) = rho * udir;
    fluxes(i, j, k, idx_t::UMX)  = rho * ux * udir + P * vdir[0];
    fluxes(i, j, k, idx_t::UMY)  = rho * uy * udir + P * vdir[1];
    fluxes(i, j, k, idx_t::UMZ)  = rho * uz * udir + P * vdir[2];
    fluxes(i, j, k, idx_t::UET)  = (rhoet + P) * udir;
  };

  // TODO: remove ParallelFor from here. Keep closures local
  void inline cons2prims(const MFIter& mfi, const Array4<Real>& cons,
                         const Array4<Real>& prims) const {

    const Box& bxg = mfi.growntilebox(idx_t::NGHOST);
    amrex::ParallelFor(bxg, [=, *this] AMREX_GPU_DEVICE(int i, int j, int k) {
      Real rho = cons(i, j, k, idx_t::URHO);
      rho = max(smallr, rho);
      Real rhoinv = Real(1.0) / rho;
      Real ux = cons(i, j, k, idx_t::UMX) * rhoinv;
      Real uy = cons(i, j, k, idx_t::UMY) * rhoinv;
      Real uz = cons(i, j, k, idx_t::UMZ) * rhoinv;
      Real rhoke = Real(0.5) * rho * (ux * ux + uy * uy + uz * uz);
      Real rhoei = cons(i, j, k, idx_t::UET) - rhoke ;
      rhoei = max(rhoei,rho*(this->get_ei_min() )); //clip energy
      Real p = (this->gamma_m1) * rhoei;

      prims(i, j, k, idx_t::QRHO) = rho;
      prims(i, j, k, idx_t::QU) = ux;
      prims(i, j, k, idx_t::QV) = uy;
      prims(i, j, k, idx_t::QW) = uz;
      prims(i, j, k, idx_t::QPRES) = p;
      prims(i, j, k, idx_t::QT) = p / (rho * this->Rspec);
      prims(i, j, k, idx_t::QC) = std::sqrt(this->gamma * p * rhoinv);
      prims(i, j, k, idx_t::QG) = this->gamma;
      prims(i, j, k, idx_t::QEINT) = rhoei * rhoinv;
      prims(i, j, k, idx_t::QFS) = 1.0;
    });
  }

  // belows are for high-order reconstruction

  AMREX_GPU_DEVICE AMREX_FORCE_INLINE void prims2cons(
      const IntVect& iv, const Array4<const Real>& prims,
      Real cons[idx_t::NCONS]) const {
    cons[idx_t::URHO] = prims(iv, idx_t::QRHO);
    cons[idx_t::UMX] = prims(iv, idx_t::QRHO) * prims(iv, idx_t::QU);
    cons[idx_t::UMY] = prims(iv, idx_t::QRHO) * prims(iv, idx_t::QV);
    cons[idx_t::UMZ] = prims(iv, idx_t::QRHO) * prims(iv, idx_t::QW);
    const Real E = max(prims(iv, idx_t::QEINT),get_ei_min() ) +
                   Real(0.5) * (prims(iv, idx_t::QU) * prims(iv, idx_t::QU) +
                                prims(iv, idx_t::QV) * prims(iv, idx_t::QV) +
                                prims(iv, idx_t::QW) * prims(iv, idx_t::QW));
    cons[idx_t::UET] = prims(iv, idx_t::QRHO) * E;
  }

  AMREX_GPU_DEVICE AMREX_FORCE_INLINE void prims2flux(
      const IntVect& iv, const int dir, const Array4<const Real>& prims,
      Real f[idx_t::NCONS]) const {
    const int QUN = idx_t::QU + dir;

    f[idx_t::URHO] = prims(iv, idx_t::QRHO) * prims(iv, QUN);
    f[idx_t::UMX] =
        prims(iv, idx_t::QRHO) * prims(iv, idx_t::QU) * prims(iv, QUN) +
        (dir == 0 ? prims(iv, idx_t::QPRES) : 0.0);
    f[idx_t::UMY] =
        prims(iv, idx_t::QRHO) * prims(iv, idx_t::QV) * prims(iv, QUN) +
        (dir == 1 ? prims(iv, idx_t::QPRES) : 0.0);
    f[idx_t::UMZ] =
        prims(iv, idx_t::QRHO) * prims(iv, idx_t::QW) * prims(iv, QUN) +
        (dir == 2 ? prims(iv, idx_t::QPRES) : 0.0);
    const Real E = prims(iv, idx_t::QEINT) +
                   Real(0.5) * (prims(iv, idx_t::QU) * prims(iv, idx_t::QU) +
                                prims(iv, idx_t::QV) * prims(iv, idx_t::QV) +
                                prims(iv, idx_t::QW) * prims(iv, idx_t::QW));
    f[idx_t::UET] =
        (prims(iv, idx_t::QPRES) + prims(iv, idx_t::QRHO) * E) * prims(iv, QUN);
  }

  /// @brief Compute the maximum characteristic speed within a stencil for local
  /// Lax-Friedrichs splitting.
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE Real
  max_char_speed(const IntVect& iv, const int dir, const int ng,
                 const Array4<const Real>& prims) const {
    const int QUN = idx_t::QU + dir;
    const amrex::IntVect ivd(amrex::IntVect::TheDimensionVector(dir));
    Real alpha = 0.0;
    for (int m = -ng; m < ng; ++m) {
      alpha = std::max(alpha, std::abs(prims(iv + m * ivd, QUN)) +
                                  prims(iv + m * ivd, idx_t::QC));
    }
    return alpha;
  }

  /// @brief Roe-averaged states between i-1 and i.
  AMREX_GPU_DEVICE
  struct RoeAvgState {
    // indices
    const int C1 = 0, C2 = 1, C3 = 2, C4 = 3;
    int CN, CT, CTT;
    // interface states
    amrex::Real u, v, w, q2, H, h, c;
  };

  /// @brief Compute the Roe-averaged states.
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE RoeAvgState roe_avg_state(
      const IntVect& iv, const int dir, const Array4<const Real>& prims) const {
    const amrex::IntVect ivm = iv - amrex::IntVect::TheDimensionVector(dir);
    RoeAvgState r;

    r.CN = dir;
    r.CT = (dir + 1) % 3;
    r.CTT = (dir + 2) % 3;

    const Real rl = prims(ivm, idx_t::QRHO);
    const Real rr = prims(iv, idx_t::QRHO);
    const Real rratio = std::sqrt(rl) / (std::sqrt(rl) + std::sqrt(rr));
    r.u = prims(ivm, idx_t::QU + r.CN) * rratio +
          prims(iv, idx_t::QU + r.CN) * (1.0 - rratio);
    r.v = prims(ivm, idx_t::QU + r.CT) * rratio +
          prims(iv, idx_t::QU + r.CT) * (1.0 - rratio);
    r.w = prims(ivm, idx_t::QU + r.CTT) * rratio +
          prims(iv, idx_t::QU + r.CTT) * (1.0 - rratio);
    r.q2 = r.u * r.u + r.v * r.v + r.w * r.w;
    const Real El = prims(ivm, idx_t::QEINT) +
                    Real(0.5) * (prims(ivm, idx_t::QU) * prims(ivm, idx_t::QU) +
                                 prims(ivm, idx_t::QV) * prims(ivm, idx_t::QV) +
                                 prims(ivm, idx_t::QW) * prims(ivm, idx_t::QW));
    const Real Er = prims(iv, idx_t::QEINT) +
                    Real(0.5) * (prims(iv, idx_t::QU) * prims(iv, idx_t::QU) +
                                 prims(iv, idx_t::QV) * prims(iv, idx_t::QV) +
                                 prims(iv, idx_t::QW) * prims(iv, idx_t::QW));
    r.H = (El + prims(ivm, idx_t::QPRES) / rl) * rratio +
          (Er + prims(iv, idx_t::QPRES) / rr) * (1.0 - rratio);
    r.h = r.H - 0.5 * r.q2;
    r.c = std::sqrt((this->gamma - 1) * r.h);

    return r;
  }

  AMREX_GPU_DEVICE AMREX_FORCE_INLINE void cons2char(
      RoeAvgState r, Real f[idx_t::NCONS]) const {
    amrex::Real tmp[idx_t::NCONS];
    amrex::Real invh = amrex::Real(1.0) / r.h;

    const int UN  = idx_t::UMX + r.CN;
    const int UT  = idx_t::UMX + r.CT;
    const int UTT = idx_t::UMX + r.CTT;

    tmp[0] = f[UT] - r.v * f[idx_t::URHO];
    tmp[1] = f[UTT] - r.w * f[idx_t::URHO];
    tmp[2] = (0.5 * ((-r.h - r.c * r.u) / r.c * f[UN] - r.v * f[UT] -
                     r.w * f[UTT] + f[idx_t::UET]) +
              (2 * r.h * r.u + r.c * r.q2) / (4 * r.c) * f[idx_t::URHO]) *
             invh;
    tmp[3] = (0.5 * ((r.h - r.c * r.u) / r.c * f[UN] - r.v * f[UT] -
                     r.w * f[UTT] + f[idx_t::UET]) +
              (-2 * r.h * r.u + r.c * r.q2) / (4 * r.c) * f[idx_t::URHO]) *
             invh;
    tmp[4] = (r.u * f[UN] + r.v * f[UT] + r.w * f[UTT] - f[idx_t::UET] +
              (r.h - 0.5 * r.q2) * f[idx_t::URHO]) *
             invh;

    for (int n = 0; n < idx_t::NCONS; ++n) f[n] = tmp[n];
  }

  AMREX_GPU_DEVICE AMREX_FORCE_INLINE void char2cons(
      RoeAvgState r, Real f[idx_t::NCONS]) const {
    amrex::Real tmp[idx_t::NCONS];

    const int UN = idx_t::UMX + r.CN;
    const int UT = idx_t::UMX + r.CT;
    const int UTT = idx_t::UMX + r.CTT;

    tmp[UN] = f[2] * (r.u - r.c) + f[3] * (r.u + r.c) + r.u * f[4];
    tmp[UT] = f[0] + r.v * (f[2] + f[3]) + r.v * f[4];
    tmp[UTT] = f[1] + r.w * (f[2] + f[3]) + r.w * f[4];
    tmp[idx_t::UET] = r.v * f[0] + r.w * f[1] + (r.H - r.u * r.c) * f[2] +
                      (r.H + r.u * r.c) * f[3] + 0.5 * r.q2 * f[4];
    tmp[idx_t::URHO] = f[2] + f[3] + f[4];

    for (int n = 0; n < idx_t::NCONS; ++n) f[n] = tmp[n];
  }

  /*
  * @brief compute primitive from conservative (local)
  * @param array of conservative vars
  * @param array of primitive variables (pass-by-pointer)
  */
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE void cons2prims_point(
    Real U[idx_t::NCONS], Real* Q) const {

    Real rho = U[idx_t::URHO];Real one_over_rho = 1.0/rho;
    Q[idx_t::QRHO] = rho;
    Real ux   = one_over_rho*U[idx_t::UMX];
    Real uy   = one_over_rho*U[idx_t::UMY];
    Real uz   = one_over_rho*U[idx_t::UMZ];
    // T and P
    Real rhoke = Real(0.5) * rho * (ux * ux+ uy* uy + uz * uz);
    Real rhoei = U[idx_t::UET] - rhoke;
    rhoei = max(rhoei,rho*(this->get_ei_min())); //clip energy
    Real p = (this->gamma_m1) * rhoei;
    Real T = p / (rho * this->Rspec);
    Q[idx_t::QT] = T;
    Q[idx_t::QPRES] = p;
    Q[idx_t::QU] = ux;
    Q[idx_t::QV] = uy;
    Q[idx_t::QW] = uz;
    // aux primitives
    Q[idx_t::QC] =  std::sqrt(gamma * p*one_over_rho);
    Q[idx_t::QG] = gamma; 
    Q[idx_t::QEINT] = cv*T;
  }

};

#endif
