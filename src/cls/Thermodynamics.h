#ifndef THERMODYNAMICS_H_
#define THERMODYNAMICS_H_

using namespace amrex;

#include <Constants.h>
#include <CNSconstants.h>


using namespace universal_constants;
using namespace CNSConstants;

// Philosophy
// Keep all thermodynamics functions local, acting on a single point

////////////////////////////////////////////////////////////////////////////////////////
template <typename idx_t>
class calorifically_perfect_gas_t {
 protected:
 public:
  Real gamma   = 1.40;   // ratio of specific heats
  Real mw = 28.96e-3;  // mean molecular weight air kg/mol

  Real gamma_m1 = gamma - Real(1.0);
  Real Ru = gas_constant;
  Real cv = Ru / (mw * gamma_m1);
  Real cp = gamma * cv;
  Real Rspec = Ru / mw;
#if CLIP_TEMPERATURE_MIN  
  Real ei_min = Rspec*min_euler_temp/gamma_m1;    // if physical T used
#else  
  Real ei_min = 1.0*min_euler_press/gamma_m1;     // if no physical T used
#endif  

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void RYP2E(const Real R,
                                                      const Real* /*Y*/,
                                                      const Real P,
                                                      Real& E) const {
    E = P / (R * gamma_m1);
  }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void RYE2TP(const Real R,
                                                       const Real* /*Y*/,
                                                       const Real E, Real& T,
                                                       Real& P) const {
    P = gamma_m1 * R * E;
    T = P / (R * Rspec);
  }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void RYE2Cs(const Real /*R*/,
                                                       const Real* /*Y*/,
                                                       const Real E,
                                                       Real& cs) const {
    Real T = E / cv;
    cs = std::sqrt(gamma * Rspec * T);
  }

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
    eint = max(eint,ei_min); //clip energy
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
      rhoei = max(rhoei,rho*(this->ei_min)); //clip energy
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
    const Real E = max(prims(iv, idx_t::QEINT),ei_min) +
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
    rhoei = max(rhoei,rho*(this->ei_min)); //clip energy
    Real p = (this->gamma_m1) * rhoei;
    Q[idx_t::QT] = p / (rho * this->Rspec);
    Q[idx_t::QPRES] = p;
  }

};


  

////////////////////////////////////////////////////////////////////////////////////////
// TODO 
class calorifically_perfect_gas_nasg_liquid_t {
 private:
  /* data */
  // gamma_a =;
  // gamma_l =;

 public:
  // state

  Real inline energy(Real p, Real rho) {
    Real eint = 0;
    eint = p/rho; // temp just to avoid warnings
    return eint;
  }

  // pressure (eint,rho)
};

#ifdef USE_PELEPHYSICS
#include <PelePhysics.H>

// Wrapper for PelePhysics EoS. No other places of the code should call pele::physics::PhysicsType::eos().
////////////////////////////////////////////////////////////////////////////////////////
template <typename idx_t>
class multispecies_pele_gas_t {
 protected:
  idx_t idx;

 public:
  
  /**
  *  @brief compute internal energy from density, pressure and mass fraction
  *  @param R density  (kg/m3)
  *  @param P pressure (Pa)
  *  @param Y array of mass fractions of chemical species
  */
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void RYP2E(const Real R,
                                                      const Real Y[NUM_SPECIES],
                                                      const Real P,
                                                      Real& E) const {
    // SI -> CGS
    Real rho_cgs = R * rho_si2cgs;
    Real p_cgs = P * pres_si2cgs;

    auto eos = pele::physics::PhysicsType::eos();
    Real e_cgs;
    eos.RYP2E(rho_cgs, Y, p_cgs, e_cgs);
    // CGS -> SI
    E = e_cgs * specenergy_cgs2si;
  }

  /**
  *  @brief compute temperature and pressure from density, internal energy and mass fraction
  *  @param R density (kg/m3)
  *  @param E mass specific internal energy  (J/kg)
  *  @param Y array of mass fractions of chemical species
  *  @return P and T
  */
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void RYE2TP(
      const Real R, const Real Y[NUM_SPECIES], const Real E, Real& T,
      Real& P) const {
    // SI -> CGS
    Real rho_cgs = R * rho_si2cgs;
    Real e_cgs = E * specenergy_si2cgs;

    auto eos = pele::physics::PhysicsType::eos();
    Real p_cgs;
    T = 0.0;
    eos.REY2T(rho_cgs, e_cgs, Y, T);
    eos.RTY2P(rho_cgs, T, Y, p_cgs);

    // CGS -> SI
    P = p_cgs * pres_cgs2si;

#ifdef DEBUG_PP_EOS
    AMREX_ALWAYS_ASSERT(T > 0.0);
    AMREX_ALWAYS_ASSERT(P > 0.0);
#endif
  }

  /**
  *  @brief compute speed of sound from density, internal energy and mass fraction
  *  @param R density (kg/m3)
  *  @param E mass specific internal energy  (J/kg)
  *  @param Y array of mass fractions of chemical species
  */
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void RYE2Cs(
      const Real R, const Real Y[NUM_SPECIES], const Real E, Real& cs) const {
    // SI -> CGS
    Real rho_cgs = R * rho_si2cgs;
    Real e_cgs = E * specenergy_si2cgs;

    auto eos = pele::physics::PhysicsType::eos();
    Real T = 0.0, p_cgs, G, cs_cgs;
    eos.REY2T(rho_cgs, e_cgs, Y, T);
    eos.RTY2P(rho_cgs, T, Y, p_cgs);
    eos.RTY2G(rho_cgs, T, Y, G);
    cs_cgs = std::sqrt(G * p_cgs / rho_cgs);
    // CGS -> SI
    cs = cs_cgs * speed_cgs2si;
#ifdef DEBUG_PP_EOS
    AMREX_ALWAYS_ASSERT(cs > 0.0);
#endif
  }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void RYE2TPCsG(
      const Real R, const Real Y[NUM_SPECIES], const Real E, Real& T, Real& P,
      Real& cs, Real& G) const {
#ifdef DEBUG_PP_EOS
    AMREX_ALWAYS_ASSERT(R > 0.0);
    AMREX_ALWAYS_ASSERT(Y[0] >= 0.0);
#endif
    // SI -> CGS
    Real rho_cgs = R * rho_si2cgs;
    Real e_cgs = E * specenergy_si2cgs;

    auto eos = pele::physics::PhysicsType::eos();
    Real p_cgs, cs_cgs;
    T = 0.0;
    eos.REY2T(rho_cgs, e_cgs, Y, T);
    eos.RTY2P(rho_cgs, T, Y, p_cgs);
    eos.RTY2G(rho_cgs, T, Y, G);
    cs_cgs = std::sqrt(G * p_cgs / rho_cgs);

    // CGS -> SI
    P = p_cgs * pres_cgs2si;
    cs = cs_cgs * speed_cgs2si;
#ifdef DEBUG_PP_EOS
    AMREX_ALWAYS_ASSERT(P > 0.0);
    AMREX_ALWAYS_ASSERT(cs > 0.0);
#endif
  }

  /**
  *  @brief compute density from pressure,temperature and mass fractions
  *  @param P pressure (kg/m3)
  *  @param T Temperature  (K)
  *  @param Y array of mass fractions of chemical species
  */
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void PYT2R(
      const Real P, const Real Y[NUM_SPECIES], const Real T, Real& R) const {
 
      auto eos = pele::physics::PhysicsType::eos();                   
      eos.PYT2R(P*pres_si2cgs, Y, T, R); R = R*rho_cgs2si;
  }            

  //-------------------------------------------------------------------------------------
  // @brief Compute mole fraction array from mass fraction species array
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void Y2X(const Real Y[NUM_SPECIES], Real (&X)[NUM_SPECIES]){
    auto eos = pele::physics::PhysicsType::eos();
    eos.Y2X(Y,X);
  }
  // @brief Compute enthalpy species array from mass fraction species array, T and rho
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void RTY2Hi(const Real rho, const Real T,const Real Y[NUM_SPECIES],
    Real (&hk)[NUM_SPECIES]) {
    auto eos = pele::physics::PhysicsType::eos();
    eos.RTY2Hi(rho, T, Y, hk);
  } 
  //-------------------------------------------------------------------------------------
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE GpuArray<Real, idx_t::NWAVES>
  cons2eigenvals(const int i, const int j, const int k,
                 const Array4<Real>& cons, const GpuArray<int, 3>& vdir) const {
   // auto eos = pele::physics::PhysicsType::eos();

    // 1. Compute density
    Real rho = 0.0;
    for (int n = 0; n < NUM_SPECIES; ++n) {
      rho += cons(i, j, k, idx_t::UFS + n);
    }

    // 2. Compute sound speed
    // Get composition (Y), energy (ei), temperature (T), pressure (p), gamma (G)
    Real rhoinv = Real(1.0) / rho;
    Real Y[NUM_SPECIES];
    for (int n = 0; n < NUM_SPECIES; ++n) {
      Y[n] = cons(i, j, k, idx_t::UFS + n) * rhoinv;
    }
    GpuArray<Real, AMREX_SPACEDIM> vel = {AMREX_D_DECL(
        cons(i, j, k, idx_t::UMX) * rhoinv, cons(i, j, k, idx_t::UMY) * rhoinv,
        cons(i, j, k, idx_t::UMZ) * rhoinv)};
    Real ke = Real(0.5) * rho *
              AMREX_D_TERM(vel[0] * vel[0], +vel[1] * vel[1], +vel[2] * vel[2]);
    Real ei = (cons(i, j, k, idx_t::UET) - ke) / rho;
    Real T, p, cs, gamma;
    this->RYE2TPCsG(rho, Y, ei, T, p, cs, gamma);

    // 3. Directional velocity abs(u)
    Real u =
        AMREX_D_TERM(vel[0] * vdir[0], +vel[1] * vdir[1], +vel[2] * vdir[2]);
    u = std::abs(u);
    // 4. Compute eigenvalues
    GpuArray<Real, idx_t::NWAVES> eigenvals;  //{u+cs,u...,u,u-cs};
    eigenvals[0] = u + cs;
    for (int n = 1; n < idx_t::NWAVES - 1; ++n) {
      eigenvals[n] = u;
    }
    eigenvals[idx_t::NWAVES-1] = u + cs;
    return eigenvals;
  }


  // void prims2char(){};

  AMREX_GPU_DEVICE AMREX_FORCE_INLINE void prims2fluxes(
      int& i, int& j, int& k, const Array4<Real>& prims, Array4<Real>&
      fluxes, const GpuArray<int, 3>& vdir) {

    Real rho = prims(i, j, k, idx.QRHO);
    Real ux = prims(i, j, k, idx.QU);
    Real uy = prims(i, j, k, idx.QV);
    Real uz = prims(i, j, k, idx.QW);
    Real P = prims(i, j, k, idx.QPRES);
    Real udir = ux * vdir[0] + uy * vdir[1] + uz * vdir[2];

    Real ekin = Real(0.5) * (ux * ux + uy * uy + uz * uz);
    Real ei;
    this->RYP2E(rho, &prims(i, j, k, idx.QFS), P, ei);
    Real rhoet = rho * ei + ekin;

    fluxes(i, j, k, idx.URHO) = rho * udir;
    fluxes(i, j, k, idx.UMX) = rho * ux * udir + P * vdir[0];
    fluxes(i, j, k, idx.UMY) = rho * uy * udir + P * vdir[1];
    fluxes(i, j, k, idx.UMZ) = rho * uz * udir + P * vdir[2];
    fluxes(i, j, k, idx.UET) = (rhoet + P) * udir;
    for (int n = 0; n < NUM_SPECIES; ++n) {
      fluxes(i, j, k, idx.UFS + n) = rho * prims(i, j, k, idx.QFS + n) * udir;
    }
  };

  // can move this to closures derived tyoe (closures_dt)
  // prims to cons
  // - We want to call it from thermodynamics class
  // - cls is stored on cpu and gpu
  // TODO: remove ParallelFor from here. Keep closures local
  void inline cons2prims(const MFIter& mfi, const Array4<Real>& cons,
                         const Array4<Real>& prims) const {
    const Box& bxg = mfi.growntilebox(idx_t::NGHOST);

    amrex::ParallelFor(bxg, [=, *this] AMREX_GPU_DEVICE(int i, int j, int k) {
      Real rho = 0.0;
      for (int n = 0; n < NUM_SPECIES; ++n) {
        rho += cons(i, j, k, idx.UFS + n);
      }
      Real rhoinv = Real(1.0) / rho;
      Real ux = cons(i, j, k, idx.UMX) * rhoinv;
      Real uy = cons(i, j, k, idx.UMY) * rhoinv;
      Real uz = cons(i, j, k, idx.UMZ) * rhoinv;
      Real rhoke = Real(0.5) * rho * (ux * ux + uy * uy + uz * uz);
      Real ei = (cons(i, j, k, idx.UET) - rhoke) * rhoinv;
      Real Y[NUM_SPECIES];
      for (int n = 0; n < NUM_SPECIES; ++n) {
        Y[n] = cons(i, j, k, idx.UFS + n) * rhoinv;
      }
      Real T, p, cs, gamma;
      this->RYE2TPCsG(rho, Y, ei, T, p, cs, gamma);

      prims(i, j, k, idx.QRHO) = rho;
      prims(i, j, k, idx.QU) = ux;
      prims(i, j, k, idx.QV) = uy;
      prims(i, j, k, idx.QW) = uz;
      prims(i, j, k, idx.QPRES) = p;
      prims(i, j, k, idx.QT) = T;
      prims(i, j, k, idx.QC) = cs;
      prims(i, j, k, idx.QG) = gamma; // why this is needed?
      prims(i, j, k, idx.QEINT) = ei;
      for (int n = 0; n < NUM_SPECIES; ++n) {
        prims(i, j, k, idx.QFS + n) = Y[n];
      }
#ifdef DEBUG_PP_EOS
      AMREX_ALWAYS_ASSERT(prims(i, j, k, idx.QPRES) > 0.0 && isnan(prims(i, j, k, idx.QPRES)) == 0);
      AMREX_ALWAYS_ASSERT(prims(i, j, k, idx.QRHO) > 0.0 && isnan(prims(i, j, k, idx.QRHO)) == 0);
      AMREX_ALWAYS_ASSERT(prims(i, j, k, idx.QT) > 0.0 && !isnan(prims(i, j, k, idx.QT)));
      AMREX_ALWAYS_ASSERT(prims(i, j, k, idx.QC) > 0.0 && isnan(prims(i, j, k, idx.QC)) == 0);
      AMREX_ALWAYS_ASSERT(prims(i, j, k, idx.QG) > 0.0 && isnan(prims(i, j, k, idx.QG)) == 0);
      AMREX_ALWAYS_ASSERT(prims(i, j, k, idx.QFS) >= 0.0 && isnan(prims(i, j, k, idx.QFS)) == 0);
      AMREX_ALWAYS_ASSERT(prims(i, j, k, idx.QFS + NUM_SPECIES - 1) >= 0.0 && isnan(prims(i, j, k, idx.QFS + NUM_SPECIES - 1)) == 0);
#endif
    });
  }

  /*
  * @brief compute primitive from conservative (local)
  * @param array of conservative vars   [NCONS]
  * @param array of primitive variables [NPRIM]
  */
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE void cons2prims_point(
    Real U[idx_t::NCONS], Real* Q ) {

    Real rho=0.0;
    for (int n = 0; n < NUM_SPECIES; ++n){
      rho += U[idx.UFS + n ];
    }
    Real rhoinv = Real(1.0) / rho;    
    Real ux = U[idx.UMX] * rhoinv;
    Real uy = U[idx.UMY] * rhoinv;
    Real uz = U[idx.UMZ] * rhoinv;
    Real rhoke = Real(0.5) * rho * (ux * ux + uy * uy + uz * uz);
    Real ei = (U[idx.UET] - rhoke) * rhoinv;
    Real Y[NUM_SPECIES];
    for (int n = 0; n < NUM_SPECIES; ++n) {
      Y[n] = U[idx.UFS + n] * rhoinv;
    }
    Real T, p, cs, gamma;
    this->RYE2TPCsG(rho, Y, ei, T, p, cs, gamma);
    Q[idx.QRHO] = rho;
    Q[idx.QU] = ux;
    Q[idx.QV] = uy;
    Q[idx.QW] = uz;
    Q[idx.QPRES] = p;
    Q[idx.QT] = T;
    Q[idx.QC] = cs;
    Q[idx.QG] = gamma; // why this is needed?
    Q[idx.QEINT] = ei;
    for (int n = 0; n < NUM_SPECIES; ++n) {
      Q[idx.QFS + n] = Y[n];
    }  
  }
  

  // belows are for high-order reconstruction

  AMREX_GPU_DEVICE AMREX_FORCE_INLINE void prims2cons(
      const IntVect& iv, const Array4<const Real>& prims,
      Real cons[idx_t::NCONS]) const {
    cons[idx_t::UMX] = prims(iv, idx_t::QRHO) * prims(iv, idx_t::QU);
    cons[idx_t::UMY] = prims(iv, idx_t::QRHO) * prims(iv, idx_t::QV);
    cons[idx_t::UMZ] = prims(iv, idx_t::QRHO) * prims(iv, idx_t::QW);
    const Real E = prims(iv, idx_t::QEINT) +
                   Real(0.5) * (prims(iv, idx_t::QU) * prims(iv, idx_t::QU) +
                                prims(iv, idx_t::QV) * prims(iv, idx_t::QV) +
                                prims(iv, idx_t::QW) * prims(iv, idx_t::QW));
    cons[idx_t::UET] = prims(iv, idx_t::QRHO) * E;
    for (int n = 0; n < NUM_SPECIES; ++n) {
      cons[idx_t::UFS + n] = prims(iv, idx_t::QRHO) * prims(iv, idx_t::QFS + n);
    }
  }

  AMREX_GPU_DEVICE AMREX_FORCE_INLINE void prims2flux(
      const IntVect& iv, const int dir, const Array4<const Real>& prims,
      Real f[idx_t::NCONS]) const {
    const int QUN = idx_t::QU + dir;

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
    for (int n = 0; n < NUM_SPECIES; ++n) {
      f[idx_t::UFS + n] =
          prims(iv, idx_t::QRHO) * prims(iv, idx_t::QFS + n) * prims(iv, QUN);
    }
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
    amrex::Real u, v, w, q2, H, h, c, Y[NUM_SPECIES];
  };

  /// @brief Compute the Roe-averaged states.
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE RoeAvgState roe_avg_state(
      const IntVect& iv, const int dir, const Array4<const Real>& prims) const {
    const amrex::IntVect ivm(iv - amrex::IntVect::TheDimensionVector(dir));
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
    r.c =
        prims(ivm, idx_t::QC) * rratio + prims(iv, idx_t::QC) * (1.0 - rratio);
    for (int ns = 0; ns < NUM_SPECIES; ++ns) {
      r.Y[ns] = prims(ivm, idx_t::QFS + ns) * rratio +
                prims(iv, idx_t::QFS + ns) * (1.0 - rratio);
    }

    return r;
  }

  AMREX_GPU_DEVICE AMREX_FORCE_INLINE void cons2char(
      RoeAvgState r, Real f[idx_t::NCONS]) const {
    amrex::Real tmp[idx_t::NCONS];
    amrex::Real invh = amrex::Real(1.0) / r.h;

    const int UN = idx_t::UMX + r.CN;
    const int UT = idx_t::UMX + r.CT;
    const int UTT = idx_t::UMX + r.CTT;

    tmp[0] = f[UT];
    for (int ns = 0; ns < NUM_SPECIES; ++ns) {
      tmp[0] += -r.v * f[idx_t::UFS + ns];
    }

    tmp[1] = f[UTT];
    for (int ns = 0; ns < NUM_SPECIES; ++ns) {
      tmp[1] += -r.w * f[idx_t::UFS + ns];
    }

    tmp[2] = 0.5 * ((-r.h - r.c * r.u) / r.c * f[UN] - r.v * f[UT] -
                    r.w * f[UTT] + f[idx_t::UET]);
    for (int ns = 0; ns < NUM_SPECIES; ++ns) {
      tmp[2] += (2 * r.h * r.u + r.c * r.q2) / (4 * r.c) * f[idx_t::UFS + ns];
    }
    tmp[2] *= invh;

    tmp[3] = 0.5 * ((r.h - r.c * r.u) / r.c * f[UN] - r.v * f[UT] -
                    r.w * f[UTT] + f[idx_t::UET]);
    for (int ns = 0; ns < NUM_SPECIES; ++ns) {
      tmp[3] += (-2 * r.h * r.u + r.c * r.q2) / (4 * r.c) * f[idx_t::UFS + ns];
    }
    tmp[3] *= invh;

    for (int ns = 0; ns < NUM_SPECIES; ++ns) {
      tmp[4 + ns] = r.Y[ns] * (r.u * f[UN] + r.v * f[UT] + r.w * f[UTT] -
                               f[idx_t::UET]) +
                    r.h * f[idx_t::UFS + ns];
      for (int n1 = 0; n1 < NUM_SPECIES; ++n1) {
        tmp[4 + ns] += -r.Y[ns] * 0.5 * r.q2 * f[idx_t::UFS + n1];
      }
      tmp[4 + ns] *= invh;
    }

    for (int n = 0; n < idx_t::NCONS; ++n) f[n] = tmp[n];
  }

  AMREX_GPU_DEVICE AMREX_FORCE_INLINE void char2cons(
      RoeAvgState r, Real f[idx_t::NCONS]) const {
    amrex::Real tmp[idx_t::NCONS];

    const int UN = idx_t::UMX + r.CN;
    const int UT = idx_t::UMX + r.CT;
    const int UTT = idx_t::UMX + r.CTT;

    tmp[UN] = f[2] * (r.u - r.c) + f[3] * (r.u + r.c);
    for (int ns = 0; ns < NUM_SPECIES; ++ns) {
      tmp[UN] += r.u * f[4 + ns];
    }

    tmp[UT] = f[0] + r.v * (f[2] + f[3]);
    for (int ns = 0; ns < NUM_SPECIES; ++ns) {
      tmp[UT] += r.v * f[4 + ns];
    }

    tmp[UTT] = f[1] + r.w * (f[2] + f[3]);
    for (int ns = 0; ns < NUM_SPECIES; ++ns) {
      tmp[UTT] += r.w * f[4 + ns];
    }

    tmp[idx_t::UET] = r.v * f[0] + r.w * f[1] + (r.H - r.u * r.c) * f[2] +
                      (r.H + r.u * r.c) * f[3];
    for (int ns = 0; ns < NUM_SPECIES; ++ns) {
      tmp[idx_t::UET] += 0.5 * r.q2 * f[4 + ns];
    }

    for (int ns = 0; ns < NUM_SPECIES; ++ns) {
      tmp[idx_t::UFS + ns] = (f[2] + f[3]) * r.Y[ns] + f[4 + ns];
    }

    for (int n = 0; n < idx_t::NCONS; ++n) f[n] = tmp[n];
  }
};
#endif

#endif
