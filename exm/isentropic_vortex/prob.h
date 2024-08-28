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
typedef closures_dt<indicies_t, visc_suth_t, cond_suth_t,
                    calorifically_perfect_gas_t<indicies_t>>
    ProbClosures;
typedef rhs_dt<keep_euler_t<false, false, 4, ProbClosures>, no_diffusive_t,
               no_source_t>
    ProbRHS;

// Numerical operators
void inline inputs() {
  ParmParse pp;

  pp.add("cns.order_rk", 3);  // -2, 1, 2 or 3"
  pp.add("cns.stages_rk", 4); // 1, 2 or 3
}
////////////////////////////////////////////////////////////////////////////////

// problem parameters
struct ProbParm {
  amrex::Real mach = 0.05;
  amrex::Real beta = 0.02; // vortex strength
  amrex::Real p0 = 1e6;    // [erg cm^-3]
  amrex::Real T0 = 300.0;  // [K]
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

  // Vortex functions (VI1)
  const Real xc = 5.0;
  const Real yc = 5.0; // vortex initial pos
  const Real R = 0.5;  // vortex radius
  const Real rsq = (x - xc) * (x - xc) + (y - yc) * (y - yc);
  
  Real u0 = pparm.mach * sqrt(cls.cp * pparm.T0);

  Real u[AMREX_SPACEDIM];
  u[0] = u0 * (1 - pparm.beta * (y - yc) / R * exp(-0.5 * rsq / R / R));
  u[1] = u0 * pparm.beta * (x - xc) / R * exp(-0.5 * rsq / R / R);
  u[2] = 0.0;
  Real T = pparm.T0 -
      0.5 * (u0 * pparm.beta * u0 * pparm.beta) / cls.cp *
        exp(-rsq / R / R);

  // printf("u = %f %f %f\n", u[0], u[1], u[2]);

  Real rho  = pparm.p0 / (cls.Rspec * T);
  Real eint = cls.cv * T;

  // Set the state
  state(i, j, k, ProbClosures::URHO) = rho;
  state(i, j, k, ProbClosures::UMX)  = rho * u[0];
  state(i, j, k, ProbClosures::UMY)  = rho * u[1];
  state(i, j, k, ProbClosures::UMZ)  = rho * u[2];
  state(i, j, k, ProbClosures::UET)  =
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
