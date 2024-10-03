#ifndef THERMODYNAMICS_H_
#define THERMODYNAMICS_H_

using namespace amrex;

#include <Constants.h>

using namespace universal_constants;

// Philosophy
// Keep all thermodynamics functions local, acting on a single point

template <typename idx_t>
class calorifically_perfect_gas_t {
 protected:

 public:
  Real gamma = 1.40;   // ratio of specific heats
  Real mw = 28.96e-3;  // mean molecular weight air kg/mol
  Real Ru = gas_constant;
  Real cv = Ru / (mw * (gamma - Real(1.0)));
  Real cp = gamma * Ru / (mw * (gamma - Real(1.0)));
  Real Rspec = Ru / mw;

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void RYP2E(const Real R,
                                                      const Real* /*Y*/,
                                                      const Real P,
                                                      Real& E) const {
    E = P / (R * (gamma - Real(1.0)));
  }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void RYE2TP(const Real R,
                                                       const Real* /*Y*/,
                                                       const Real E, Real& T,
                                                       Real& P) const {
    P = (gamma - Real(1.0)) * R * E;
    T = P / (R * Rspec);
  }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void RYE2Cs(const Real /*R*/,
                                                       const Real* /*Y*/,
                                                       const Real E,
                                                       Real& cs) const {
    Real T = E / cv;
    cs = std::sqrt(gamma * Rspec * T);
  }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE GpuArray<Real,idx_t::NWAVES> cons2eigenvals(const int i, const int j, const int k, const Array4<Real>& cons, const GpuArray<int, 3>& vdir) const {

    Real rho = cons(i, j, k, idx_t::URHO);
    Real rhoinv = Real(1.0) / rho;
    GpuArray<Real, AMREX_SPACEDIM> vel= {AMREX_D_DECL(
                                  cons(i, j, k, idx_t::UMX) * rhoinv,
                                  cons(i, j, k, idx_t::UMY) * rhoinv,
                                  cons(i, j, k, idx_t::UMZ) * rhoinv)};

    Real ke = AMREX_D_PICK(
                                vel[0]*vel[0],
                                vel[0]*vel[0] + vel[1]*vel[1], 
                                vel[0]*vel[0] + vel[1]*vel[1] + vel[2]*vel[2]);
    ke = Real(0.5)*rho*ke;
    Real eint = (cons(i,j,k,idx_t::UET) - ke)/rho;
    Real T = eint / cv;

    Real cs = std::sqrt(gamma * Rspec * T);
    Real u = AMREX_D_PICK(
              vel[0]*vdir[0],
              vel[0]*vdir[0] + vel[1]*vdir[1], 
              vel[0]*vdir[0] + vel[1]*vdir[1] + vel[2]*vdir[2]);
    GpuArray<Real,idx_t::NWAVES> eigenvals={u+cs,u,u-cs};
  return eigenvals;
  }

  // Real sos(Real& T){return std::sqrt(this->gamma * this->Rspec * T);}
  // void prims2cons(i,j,k,){};

  // void prims2chars(int i, int j, int k, const Array4<Real>& prims, Array4<Real>& chars, const GpuArray<int, AMREX_SPACEDIM>& vdir){

  // };

  // void chars2prims() {}

  AMREX_GPU_DEVICE AMREX_FORCE_INLINE void prims2fluxes(
      int i, int j, int k, const Array4<Real>& prims, const Array4<Real>& fluxes, const GpuArray<int, 3>& vdir) const {
    Real rho = prims(i, j, k, idx_t::QRHO);
    Real ux = prims(i, j, k, idx_t::QU);
    Real uy = prims(i, j, k, idx_t::QV);
    Real uz = prims(i, j, k, idx_t::QW);
    Real P = prims(i, j, k, idx_t::QPRES);
    Real udir = ux * vdir[0] + uy * vdir[1] + uz * vdir[2];

    Real ekin = Real(0.5) * (ux * ux + uy * uy + uz * uz);
    Real rhoet = rho * (cp * prims(i, j, k, idx_t::QT) + ekin);

    fluxes(i, j, k, idx_t::URHO) = rho * udir;
    fluxes(i, j, k, idx_t::UMX) = rho * ux * udir + P * vdir[0];
    fluxes(i, j, k, idx_t::UMY) = rho * uy * udir + P * vdir[1];
    fluxes(i, j, k, idx_t::UMZ) = rho * uz * udir + P * vdir[2];
    fluxes(i, j, k, idx_t::UET) = (rhoet + P) * udir;
  };

  // TODO: remove ParallelFor from here. Keep closures local
  void inline cons2prims(const MFIter& mfi, const Array4<Real>& cons,
                         const Array4<Real>& prims) const {
    const Box& bxg = mfi.growntilebox(idx_t::NGHOST);

    amrex::ParallelFor(bxg, [=, *this] AMREX_GPU_DEVICE(int i, int j, int k) {
      Real rho = cons(i, j, k, idx_t::URHO);
      // Print() << "cons2prim"<< i << j << k << rho << std::endl;
      rho = max(1e-40, rho);
      Real rhoinv = Real(1.0) / rho;
      Real ux = cons(i, j, k, idx_t::UMX) * rhoinv;
      Real uy = cons(i, j, k, idx_t::UMY) * rhoinv;
      Real uz = cons(i, j, k, idx_t::UMZ) * rhoinv;
      Real rhoke = Real(0.5) * rho * (ux * ux + uy * uy + uz * uz);
      Real rhoei = (cons(i, j, k, idx_t::UET) - rhoke);
      Real p = (this->gamma - Real(1.0)) * rhoei;

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
};

class calorifically_perfect_gas_nasg_liquid_t {
 private:
  /* data */
  // gamma_a =;
  // gamma_l =;

 public:
  // state

  Real inline energy(Real p, Real rho) {
    Real eint = 0;
    return eint;
  }

  // pressure (eint,rho)
};

#ifdef USE_PELEPHYSICS
#include <PelePhysics.H>

// #define DEBUG_PP_EOS

// Wrapper for PelePhysics EoS. No other places of the code should call
// pele::physics::PhysicsType::eos().
template <typename idx_t>
class multispecies_gas_t {
 protected:
  idx_t idx;

 public:
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
    E = e_cgs*specenergy_cgs2si
  }

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

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void RYE2TPCsG(const Real R,
                                                     const Real Y[NUM_SPECIES],
                                                     const Real E, Real& T,
                                                     Real& P, Real& cs,
                                                     Real& G) const {
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

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE GpuArray<Real,idx_t::NWAVES> cons2eigenvals(const int i, const int j, const int k, const int dir,const Array4<Real>& cons) const {

    Abort("cons2eigenvals not implemented for multispecies_gas_t");

    // GpuArray<int, AMREX_SPACEDIM> vdir = {AMREX_D_DECL(int(dir == 0), int(dir == 1), int(dir == 2))};

    // Real rho = cons(i, j, k, idx.URHO);
    // Real rhoinv = Real(1.0) / rho;
    // GpuArray<int, AMREX_SPACEDIM> vel= {AMREX_D_DECL(
    //                               cons(i, j, k, idx_t::UMX) * rhoinv,
    //                               cons(i, j, k, idx_t::UMY) * rhoinv,
    //                               cons(i, j, k, idx_t::UMZ) * rhoinv)};
    // Real ke = Real(0.5) * rho * AMREX_D_PICK(
    //                             vel[0]*vel[0],
    //                             vel[0]*vel[0] + vel[1]*vel[1], 
    //                             vel[0]*vel[0] + vel[1]*vel[1] + vel[2]*vel[2]);
    // Real eint = (cons(i,j,k,idx_t::UET) - ke)/rho;
    // Real T = eint / cv;

    // Real cs = std::sqrt(gamma * Rspec * T);
    // Real u = AMREX_D_PICK(
    //           vel[0]*vdir[0],
    //           vel[0]*vdir[0] + vel[1]*vdir[1], 
    //           vel[0]*vdir[0] + vel[1]*vdir[1] + vel[2]*vdir[2]);

    GpuArray<Real,idx_t::NWAVES> eigenvals; //{u+cs,u,u-cs}
  return eigenvals;
  }

  // void prims2cons(i,j,k,){};

  // void prims2char(){};

  // unused
  // AMREX_GPU_DEVICE AMREX_FORCE_INLINE void prims2fluxes(
  //     int& i, int& j, int& k, const Array4<Real>& prims, Array4<Real>&
  //     fluxes, const GpuArray<int, 3>& vdir) {
  //   Real rho = prims(i, j, k, idx.QRHO);
  //   Real ux = prims(i, j, k, idx.QU);
  //   Real uy = prims(i, j, k, idx.QV);
  //   Real uz = prims(i, j, k, idx.QW);
  //   Real P = prims(i, j, k, idx.QPRES);
  //   Real udir = ux * vdir[0] + uy * vdir[1] + uz * vdir[2];

  //   Real ekin = Real(0.5) * (ux * ux + uy * uy + uz * uz);
  //   Real ei;
  //   this->RYP2E(rho, &prims(i, j, k, idx.QFS), P, ei);
  //   Real rhoet = rho * ei + ekin;

  //   fluxes(i, j, k, idx.URHO) = rho * udir;
  //   fluxes(i, j, k, idx.UMX) = rho * ux * udir + P * vdir[0];
  //   fluxes(i, j, k, idx.UMY) = rho * uy * udir + P * vdir[1];
  //   fluxes(i, j, k, idx.UMZ) = rho * uz * udir + P * vdir[2];
  //   fluxes(i, j, k, idx.UET) = (rhoet + P) * udir;
  //   for (int n = 0; n < NUM_SPECIES; ++n) {
  //     fluxes(i, j, k, idx.UFS + n) = rho * prims(i, j, k, idx.QFS + n) *
  //     udir;
  //   }
  // };

  // can move this to closures derived tyoe (closures_dt)
  // prims to cons
  // - We want to call it from thermodynamics class
  // - cls is stored on cpu and gpu
  // TODO: remove ParallelFor from here. Keep closures local
  void inline cons2prims(const MFIter& mfi, const Array4<Real>& cons,
                         const Array4<Real>& prims) const {
    const Box& bxg = mfi.growntilebox(idx_t::NGHOST);

    amrex::ParallelFor(bxg, [=, *this] AMREX_GPU_DEVICE(int i, int j, int k) {
      Real rho = 0.0;  // cons(i, j, k, idx.URHO);
      for (int n = 0; n < NUM_SPECIES; ++n) {
        rho += cons(i, j, k, idx.UFS + n);
      }
      // rho = max(1e-40, rho);
      Real rhoinv = Real(1.0) / rho;
      Real ux = cons(i, j, k, idx.UMX) * rhoinv;
      Real uy = cons(i, j, k, idx.UMY) * rhoinv;
      Real uz = cons(i, j, k, idx.UMZ) * rhoinv;
      Real rhoke = Real(0.5) * rho * (ux * ux + uy * uy + uz * uz);
      Real ei = (cons(i, j, k, idx.UET) - rhoke) * rhoinv;
      Real Y[NUM_SPECIES];
      for (int n = 0; n < NUM_SPECIES; ++n) {
        Y[n] = cons(i, j, k, idx.UFS + n) * rhoinv;
        // if (Y[n] < 0.0) {
        //   std::cout << "Y[" << n << "](" << i << ", " << j << ", " << k
        //             << ") = " << Y[n] << std::endl;
        //   std::cout << cons(i, j, k, idx.UFS + n) << " " << cons(i, j, k,
        //   idx.URHO) << std::endl;
        // }
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
      prims(i, j, k, idx.QG) = gamma;
      prims(i, j, k, idx.QEINT) = ei;
      for (int n = 0; n < NUM_SPECIES; ++n) {
        prims(i, j, k, idx.QFS + n) = Y[n];
      }
#ifdef DEBUG_PP_EOS
      // AMREX_DEVICE_PRINTF(static_cast<const char*>("cons2prims: %d %d %d\n"), i, j, k);
      AMREX_ALWAYS_ASSERT(prims(i, j, k, idx.QPRES) > 0.0 && isnan(prims(i, j, k, idx.QPRES)) == 0);
      AMREX_ALWAYS_ASSERT(prims(i, j, k, idx.QRHO) > 0.0 && isnan(prims(i, j, k, idx.QRHO)) == 0);
      AMREX_ALWAYS_ASSERT(prims(i, j, k, idx.QT) > 0.0 && !isnan(prims(i, j, k, idx.QT)));
      AMREX_ALWAYS_ASSERT(prims(i, j, k, idx.QC) > 0.0 && isnan(prims(i, j, k, idx.QC)) == 0);
      AMREX_ALWAYS_ASSERT(prims(i, j, k, idx.QG) > 0.0 && isnan(prims(i, j, k, idx.QG)) == 0);
      AMREX_ALWAYS_ASSERT(prims(i, j, k, idx.QFS) >= 0.0 && isnan(prims(i, j, k, idx.QFS)) == 0);
      AMREX_ALWAYS_ASSERT(prims(i, j, k, idx.QFS + NUM_SPECIES - 1) >= 0.0 && isnan(prims(i, j, k, idx.QFS + NUM_SPECIES - 1)) == 0);
      // if (i == -3) AMREX_DEVICE_PRINTF(static_cast<const char*>("T = %f: %i %i %i\n"), prims(i, j, k, idx.QT), i, j, k);
#endif
    });
  }
};
#endif

#endif
