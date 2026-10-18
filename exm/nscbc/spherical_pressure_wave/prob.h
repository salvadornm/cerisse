#ifndef CNS_PROB_H_
#define CNS_PROB_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <Closures.h>
#include <RHS.h>
#include <cmath>

using namespace amrex;
namespace PROB {
static_assert(AMREX_SPACEDIM == 3, "The spherical-wave test requires DIM=3");

// Lodato et al., JCP 227 (2008), section 4.4, Eq. (96).
// These dimensional parameters stay fixed when the computational box is doubled.
struct ProbParm {
  Real pressure = 101325.0;
  Real temperature = 300.0;
  Real amplitude = 0.001;
  Real radius = 0.00065; // 0.05 times the SMALL box's 0.013 m side
};

// Table 2 air properties; reuse CERISSE's perfect-gas thermodynamics.
struct AirEOS : calorifically_perfect_gas_t<indicies_t> {
  AirEOS() {
    Rspec = Ru / Real(28.9e-3);
    cv = Rspec / gamma_m1;
    cp = gamma * cv;
  }
};
struct AirTransport {
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  Real visc(Real T) const {
    return Real(1.827e-5) * std::pow(T / Real(291.15), Real(1.5)) *
           Real(291.15 + 120.0) / (T + Real(120.0));
  }
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  Real cond(Real T) const {
    const Real cp = Real(1.4) / Real(0.4) * universal_constants::gas_constant / Real(28.9e-3);
    return visc(T) * cp / Real(0.72);
  }
};
struct MethodParm {
  static constexpr int order = 4;
  static constexpr bool dissipation = false;
  static constexpr Real C2skew = 0.1, C4skew = 0.016;
};
struct ViscousParm { static constexpr int order = 4; };
using ProbClosures = closures_dt<indicies_t, AirEOS, AirTransport>;
using ProbRHS = rhs_dt<skew_t<MethodParm, ProbClosures>,
                      viscous_t<ViscousParm, ProbClosures>, no_source_t>;

inline void inputs() {
  amrex::Print() << "Lodato spherical pressure pulse: L=0.013 m, Rp=0.00065 m, delta=0.001\n"
                 << "Fourth-order skew-symmetric convection and centered diffusion\n";
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void prob_initdata(int i, int j, int k, Array4<Real> const& state,
                   GeometryData const& geomdata, ProbClosures const& cls,
                   ProbParm const& pp) {
  const auto* lo = geomdata.ProbLo();
  const auto* dx = geomdata.CellSize();
  const Real x = lo[0] + (Real(i)+Real(0.5))*dx[0];
  const Real y = lo[1] + (Real(j)+Real(0.5))*dx[1];
  const Real z = lo[2] + (Real(k)+Real(0.5))*dx[2];
  const Real p = pp.pressure * (Real(1.0) + pp.amplitude *
      std::exp(-(x*x+y*y+z*z)/(Real(2.0)*pp.radius*pp.radius)));
  Real rho, eint;
  cls.PYT2R(p, nullptr, pp.temperature, rho);
  cls.RYP2E(rho, nullptr, p, eint);
  state(i,j,k,cls.URHO) = rho;
  state(i,j,k,cls.UMX) = Real(0.0);
  state(i,j,k,cls.UMY) = Real(0.0);
  state(i,j,k,cls.UMZ) = Real(0.0);
  state(i,j,k,cls.UET) = rho*eint;
}

// Ordinary physical BC is extrapolation; NSCBC advances the ghost shell.
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void bcnormal(const Real*, Real, const Real* s_int, const Real*, Real* s_ext,
              int, int, Real, GeometryData const&, ProbClosures const&, ProbParm const&) {
  for (int n=0; n<ProbClosures::NCONS; ++n) s_ext[n] = s_int[n];
}
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void user_tagging(int, int, int, int, auto&, const auto&, const auto&,
                  ProbParm const&, int) {}
} // namespace PROB
#endif
