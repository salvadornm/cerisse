#ifndef THERMODYNAMICS_H_
#define THERMODYNAMICS_H_

using namespace amrex;

template <typename idx_t>
class calorifically_perfect_gas_t {
 protected:
  idx_t idx;

 public:
  Real gamma = 1.40;   // ratio of specific heats
  Real mw = 28.96e-3;  // mean molecular weight air kg/mol
  Real Ru = Real(8.31451);
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

  // Real sos(Real& T){return std::sqrt(this->gamma * this->Rspec * T);}

  // void prims2cons(i,j,k,){};

  void prims2char(){};

  // void prims2flux(int& i, int& j, int& k, const GpuArray<Real,NPRIM>& prims,
  // GpuArray<Real,NCONS>& fluxes, const GpuArray<int, 3>& vdir) {
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE void prims2fluxes(
      int& i, int& j, int& k, const Array4<Real>& prims, Array4<Real>& fluxes,
      const GpuArray<int, 3>& vdir) {
    Real rho = prims(i, j, k, idx.QRHO);
    Real ux = prims(i, j, k, idx.QU);
    Real uy = prims(i, j, k, idx.QV);
    Real uz = prims(i, j, k, idx.QW);
    Real P = prims(i, j, k, idx.QPRES);
    Real udir = ux * vdir[0] + uy * vdir[1] + uz * vdir[2];

    Real ekin = Real(0.5) * (ux * ux + uy * uy + uz * uz);
    Real rhoet = rho * (cp * prims(i, j, k, idx.QT) + ekin);

    fluxes(i, j, k, idx.URHO) = rho * udir;
    fluxes(i, j, k, idx.UMX) = rho * ux * udir + P * vdir[0];
    fluxes(i, j, k, idx.UMY) = rho * uy * udir + P * vdir[1];
    fluxes(i, j, k, idx.UMZ) = rho * uz * udir + P * vdir[2];
    fluxes(i, j, k, idx.UET) = (rhoet + P) * udir;
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
      Real rho = cons(i, j, k, idx.URHO);
      // Print() << "cons2prim"<< i << j << k << rho << std::endl;
      rho = max(1e-40, rho);
      Real rhoinv = Real(1.0) / rho;
      Real ux = cons(i, j, k, idx.UMX) * rhoinv;
      Real uy = cons(i, j, k, idx.UMY) * rhoinv;
      Real uz = cons(i, j, k, idx.UMZ) * rhoinv;
      Real rhoke = Real(0.5) * rho * (ux * ux + uy * uy + uz * uz);
      Real rhoei = (cons(i, j, k, idx.UET) - rhoke);
      Real p = (this->gamma - Real(1.0)) * rhoei;

      prims(i, j, k, idx.QRHO) = rho;
      prims(i, j, k, idx.QU) = ux;
      prims(i, j, k, idx.QV) = uy;
      prims(i, j, k, idx.QW) = uz;
      prims(i, j, k, idx.QPRES) = p;
      prims(i, j, k, idx.QT) = p / (rho * this->Rspec);
      prims(i, j, k, idx.QC) = std::sqrt(this->gamma * p * rhoinv);
      prims(i, j, k, idx.QG) = this->gamma;
      prims(i, j, k, idx.QEINT) = rhoei * rhoinv;
      prims(i, j, k, idx.QFS) = 1.0;
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
    Real rho_cgs = R * 0.001;
    Real p_cgs = P * 10.0;

    auto eos = pele::physics::PhysicsType::eos();
    Real e_cgs;
    eos.RYP2E(rho_cgs, Y, p_cgs, e_cgs);

    E = e_cgs * 1.0e-4;
    // AMREX_ALWAYS_ASSERT(E > 0.0);
  }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void RYE2TP(
      const Real R, const Real Y[NUM_SPECIES], const Real E, Real& T,
      Real& P) const {
    Real rho_cgs = R * 0.001;
    Real e_cgs = E * 1.0e4;

    auto eos = pele::physics::PhysicsType::eos();
    Real p_cgs;
    eos.REY2T(rho_cgs, e_cgs, Y, T);
    eos.RTY2P(rho_cgs, T, Y, p_cgs);

    P = p_cgs * 0.1;
    // AMREX_ALWAYS_ASSERT(P > 0.0);
  }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void RYE2Cs(
      const Real R, const Real Y[NUM_SPECIES], const Real E, Real& cs) const {
    Real rho_cgs = R * 0.001;
    Real e_cgs = E * 1.0e4;

    auto eos = pele::physics::PhysicsType::eos();
    Real T, p_cgs, G, cs_cgs;
    eos.REY2T(rho_cgs, e_cgs, Y, T);
    eos.RTY2P(rho_cgs, T, Y, p_cgs);
    eos.RTY2G(rho_cgs, T, Y, G);
    cs_cgs = std::sqrt(G * p_cgs / rho_cgs);

    cs = cs_cgs * 0.01;
    // AMREX_ALWAYS_ASSERT(cs > 0.0);
  }

  AMREX_GPU_DEVICE AMREX_FORCE_INLINE void RYE2TPCsG(const Real R,
                                                     const Real Y[NUM_SPECIES],
                                                     const Real E, Real& T,
                                                     Real& P, Real& cs,
                                                     Real& G) const {
    // AMREX_ALWAYS_ASSERT(R > 0.0);
    // AMREX_ALWAYS_ASSERT(Y[0] >= 0.0);
    // AMREX_ALWAYS_ASSERT(E > 0.0);
    Real rho_cgs = R * 0.001;
    Real e_cgs = E * 1.0e4;

    auto eos = pele::physics::PhysicsType::eos();
    Real p_cgs, cs_cgs;
    eos.REY2T(rho_cgs, e_cgs, Y, T);
    eos.RTY2P(rho_cgs, T, Y, p_cgs);
    eos.RTY2G(rho_cgs, T, Y, G);
    cs_cgs = std::sqrt(G * p_cgs / rho_cgs);

    P = p_cgs * 0.1;
    cs = cs_cgs * 0.01;
    // AMREX_ALWAYS_ASSERT(P > 0.0);
    // AMREX_ALWAYS_ASSERT(cs > 0.0);
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

      // AMREX_ALWAYS_ASSERT(prims(i, j, k, idx.QPRES) > 0.0);
      // AMREX_ALWAYS_ASSERT(prims(i, j, k, idx.QRHO) > 0.0);
      // AMREX_ALWAYS_ASSERT(prims(i, j, k, idx.QT) > 0.0);
      // AMREX_ALWAYS_ASSERT(prims(i, j, k, idx.QC) > 0.0);
      // AMREX_ALWAYS_ASSERT(prims(i, j, k, idx.QFS) >= 0.0);
    });
  }
};
#endif

#endif
