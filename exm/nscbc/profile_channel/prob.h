#ifndef CNS_PROB_H_
#define CNS_PROB_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_ParmParse.H>
#include <Closures.h>
#include <RHS.h>
#include <bc_types.h>

using namespace amrex;

namespace PROB {

struct ProbParm {
  Real gamma = 1.4;
  Real p0 = 101325.0;
  Real T0 = 300.0;
  Real Rair = 287.0;
  Real rho0 = p0 / (Rair*T0);
  Real c0 = std::sqrt(gamma*p0/rho0);
  Real Mach = 0.1;
  Real umax = Mach*c0;
  Real half_height = 1.0e-3;
  Real Y0[NUM_SPECIES] = {0.0};
};

struct methodparm_t {
  static constexpr bool dissipation = true;
  static constexpr int order = 4;
  static constexpr Real C2skew = 0.01;
  static constexpr Real C4skew = 0.016;
};

using ProbClosures = closures_dt<indicies_t, visc_suth_t, cond_suth_t,
                                 calorifically_perfect_gas_t<indicies_t>>;
using ProbRHS = rhs_dt<keep_euler_t<false, false, 4, ProbClosures>,
                       viscous_t<methodparm_t, ProbClosures>, no_source_t>;
using GlobalBC = manual_bc_t<ProbClosures>;

inline void inputs() {
  amrex::Print() << "2-D channel: profiled NSCBC inflow and pressure-relaxed "
                    "non-reflecting outflow\n";
}

#ifdef USE_MANUAL_NSBC_TARGET
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void nscbc_target(amrex::Real /*x*/, amrex::Real y, amrex::Real /*z*/,
                  int dir, int side_sign, amrex::Real& u,
                  amrex::Real& v, amrex::Real& w, amrex::Real& T)
{
  // x-low inlet. Coordinates are physical face coordinates.  This profile
  // is zero at y = +/-h and reaches umax on the channel centreline.
  if (dir == 0 && side_sign > 0) {
    constexpr amrex::Real h = 1.0e-3;
    constexpr amrex::Real T0 = 300.0;
    constexpr amrex::Real gamma = 1.4;
    constexpr amrex::Real Rair = 287.0;
    constexpr amrex::Real Mach = 0.1;
    const amrex::Real c0 = sqrt(gamma*Rair*T0);
    const amrex::Real yn = y/h;
    u = Mach*c0*amrex::max(amrex::Real(0.0), amrex::Real(1.0) - yn*yn);
    v = amrex::Real(0.0);
    w = amrex::Real(0.0);
    T = T0;
  }
}
#endif

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void prob_initdata(int i, int j, int k, Array4<Real> const& state,
                   GeometryData const& geomdata, ProbClosures const& cls,
                   ProbParm const& pp)
{
  const Real* lo = geomdata.ProbLo();
  const Real* dx = geomdata.CellSize();
  const Real y = lo[1] + (Real(j) + Real(0.5))*dx[1];
  const Real yn = y/pp.half_height;
  const Real u = pp.umax*amrex::max(Real(0.0), Real(1.0) - yn*yn);
  const Real rho = pp.p0/(pp.Rair*pp.T0);
  const Real rhoe = pp.p0/(cls.gamma - Real(1.0));

  state(i,j,k,cls.URHO) = rho;
  state(i,j,k,cls.UMX) = rho*u;
  state(i,j,k,cls.UMY) = Real(0.0);
  state(i,j,k,cls.UMZ) = Real(0.0);
  state(i,j,k,cls.UET) = rhoe + Real(0.5)*rho*u*u;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void bcnormal(const Real[AMREX_SPACEDIM], Real,
              const Real[ProbClosures::NCONS],
              const Real[ProbClosures::NCONS],
              Real[ProbClosures::NCONS], const int, const int,
              const Real, GeometryData const&, ProbClosures const&,
              ProbParm const&) {}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void user_tagging(int, int, int, int, auto&, const auto&, const auto&,
                  const ProbParm&, int) {}

} // namespace PROB
#endif
