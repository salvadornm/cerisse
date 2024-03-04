#ifndef PROB_H
#define PROB_H

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_ParmParse.H>
#include <Closures.h>
#include <RHS.h>

using namespace amrex;

namespace PROB {

//////////////////////////////DISCRETISATION////////////////////////////////////
inline Vector<std::string> cons_vars_names={"Density","Xmom","Ymom","Zmom","Energy"};
inline Vector<int> cons_vars_type={0,1,2,3,0};

typedef closures_dt<indicies_t, visc_suth_t, cond_suth_t,
                    calorifically_perfect_gas_t<indicies_t>>
    ProbClosures;
typedef rhs_dt<keep_euler_t<false, false, 4, ProbClosures>, no_diffusive_t,
               no_source_t>
    ProbRHS;

// Numerical operators
void inline inputs() {
  ParmParse pp;

  //-1 = N/A (Incase of periodic)
  // 0 = Interior           3 = Symmetry
  // 1 = Inflow             4 = SlipWall
  // 2 = Outflow            5 = NoSlipWall
  // 6 = user defined
  pp.addarr("cns.lo_bc", std::vector<int>{1, 1, 1});
  pp.addarr("cns.hi_bc", std::vector<int>{1, 1, 1});
  pp.add("cns.order_rk", 3);  // -2, 1, 2 or 3"
  pp.add("cns.stages_rk", 4); // 1, 2 or 3
  pp.add("cns.screen_output", 100000); //
  pp.add("cns.verbose", 0);            // 0=quiet, 1=verbose
}
////////////////////////////////////////////////////////////////////////////////

// problem parameters
struct ProbParm {
  bool convecting = false;

  Real Re = Real(1600.0);
  Real mach = Real(1.0);
  Real prandtl = Real(0.71);
  Real omega_x = Real(1.0); // [rad s^-1]
  Real omega_y = Real(1.0); // [rad s^-1]
  Real omega_z = Real(1.0); // [rad s^-1]
  Real L = Real(1.0);       // [m]
  Real T0 = Real(110.4);    // [K]
};

// initial condition
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
prob_initdata(int i, int j, int k, Array4<Real> const &state,
              GeometryData const &geomdata, ProbClosures const &cls,
              const ProbParm &pparm) {
  // Geometry
  const Real *prob_lo = geomdata.ProbLo();
  const Real *dx = geomdata.CellSize();
  const Real x = prob_lo[0] + (i + 0.5) * dx[0];
  const Real y = prob_lo[1] + (j + 0.5) * dx[1];
  const Real z = prob_lo[2] + (k + 0.5) * dx[2];

  Real cs = std::sqrt(cls.gamma * cls.Rspec * pparm.T0);
  Real v0 = pparm.mach * cs;
  // density (rho0 = Re*mu0/(L*U0)
  Real rho0 = pparm.Re * cls.visc(pparm.T0) / (pparm.L * v0);
  // pressure (p0 = rho0*R*T0)
  Real p0 = rho0 * cls.Rspec * pparm.T0;

  // TGV functions
  Real u[3] = {0.0};
  u[0] = v0 * sin(pparm.omega_x * x / pparm.L) *
         cos(pparm.omega_y * y / pparm.L) * cos(pparm.omega_z * z / pparm.L);
  u[1] = -v0 * cos(pparm.omega_x * x / pparm.L) *
         sin(pparm.omega_y * y / pparm.L) * cos(pparm.omega_z * z / pparm.L);
  if (pparm.convecting) {
    u[0] += v0;
    u[1] += v0;
  }
  const Real p = p0 + rho0 * v0 * v0 / Real(16.0) *
                          (cos(2.0 * pparm.omega_x * x / pparm.L) +
                           cos(2.0 * pparm.omega_y * y / pparm.L)) *
                          (cos(2.0 * pparm.omega_z * z / pparm.L) + Real(2.0));
  Real rho = p / (cls.Rspec * pparm.T0);
  Real eint = cls.cv * pparm.T0;

  // Set the state
  state(i, j, k, cls.URHO) = rho;
  state(i, j, k, cls.UMX) = rho * u[0];
  state(i, j, k, cls.UMY) = rho * u[1];
  state(i, j, k, cls.UMZ) = rho * u[2];
  state(i, j, k, cls.UET) =
      rho * (eint + Real(0.5) * (u[0] * u[0] + u[1] * u[1] + u[2] * u[2]));
}

// boundary conditions
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
bcnormal(const Real x[AMREX_SPACEDIM], Real dratio, const Real s_int[ProbClosures::NCONS],
         const Real s_refl[ProbClosures::NCONS], Real s_ext[ProbClosures::NCONS], const int idir,
         const int sgn, const Real time, GeometryData const & /*geomdata*/,
         ProbClosures const &cls, ProbParm const &pparm) {
  Abort("bcnormal not set");
}

// source term
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
user_source(int i, int j, int k, const auto &state, const auto &rhs,
            const ProbParm &pparm, ProbClosures const &cls, auto const dx) {}
////////////////////////////////////////////////////////////////////////////////

///////////////////////////////AMR//////////////////////////////////////////////
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
user_tagging(int i, int j, int k, int nt_level, auto &tagfab,
             const auto &sdatafab, const auto &geomdata, const ProbParm &pparm,
             int level) {}
////////////////////////////////////////////////////////////////////////////////

} // namespace PROB
#endif
