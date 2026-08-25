#ifndef CNS_PROB_H_
#define CNS_PROB_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <PelePhysics.H>
#include <ReactorBase.H>

#include <Closures.h>
#include <RHS.h>

using namespace amrex;

namespace PROB {

struct ProbParm {
  Real pressure = 101325.0;
  Real temperature = 300.0;
  Real velocity = 20.0;
  Real centre_x = 3.0e-3;
  Real centre_y = 0.0;
  Real radius = 7.5e-4;
};

struct methodparm_t {
  static constexpr bool dissipation = true;
  static constexpr int order = 4;
  static constexpr Real C2skew = 0.01;
  static constexpr Real C4skew = 0.016;
};

using ProbClosures = closures_dt< indicies_t, transport_Pele_t, multispecies_pele_gas_t<indicies_t>>;

//using ProbRHS = rhs_dt<skew_t<methodparm_t, ProbClosures>,no_diffusive_t, no_source_t>;
using ProbRHS =  rhs_dt< keep_euler_spec_t<false, false, 4, ProbClosures>, no_diffusive_t, no_source_t>;


inline void inputs() {
  amrex::Print()
      << "2-D O2/N2 Gaussian species convection with NSCBC inlet/outlet\n"
      << "PelePhysics mechanism: air; convection velocity: 20 m/s\n";
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void prob_initdata(int i, int j, int k, Array4<Real> const& state,
                   GeometryData const& geomdata, ProbClosures const& cls,
                   ProbParm const& pp) {
  const Real* lo = geomdata.ProbLo();
  const Real* dx = geomdata.CellSize();
  const Real x = lo[0] + (Real(i) + Real(0.5))*dx[0];
  const Real y = lo[1] + (Real(j) + Real(0.5))*dx[1];
  const Real r2 = (x-pp.centre_x)*(x-pp.centre_x) +
                  (y-pp.centre_y)*(y-pp.centre_y);

  Real Y[NUM_SPECIES] = {Real(0.0)};
  Y[O2_ID] = exp(-r2/(Real(2.0)*pp.radius*pp.radius));
  Y[N2_ID] = Real(1.0) - Y[O2_ID];

  Real rho = Real(0.0);
  Real eint = Real(0.0);
  cls.PYT2R(pp.pressure, Y, pp.temperature, rho);
  cls.RYP2E(rho, Y, pp.pressure, eint);

  state(i,j,k,cls.UMX) = rho*pp.velocity;
  state(i,j,k,cls.UMY) = Real(0.0);
  state(i,j,k,cls.UMZ) = Real(0.0);
  state(i,j,k,cls.UET) = rho*eint + Real(0.5)*rho*pp.velocity*pp.velocity;
  for (int ns = 0; ns < NUM_SPECIES; ++ns) {
    state(i,j,k,cls.UFS + ns) = rho*Y[ns];
  }
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void bcnormal(const Real[AMREX_SPACEDIM], Real,
              const Real[ProbClosures::NCONS],
              const Real[ProbClosures::NCONS],
              Real[ProbClosures::NCONS], int, int, Real,
              GeometryData const&, ProbClosures const&, ProbParm const&) {}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void user_tagging(int, int, int, int, auto&, const auto&, const auto&,
                  const ProbParm&, int) {}

} // namespace PROB
#endif
