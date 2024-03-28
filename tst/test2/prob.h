#ifndef CNS_PROB_H_
#define CNS_PROB_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_ParmParse.H>
#include <PelePhysics.H>
#include <ReactorBase.H>

#include "Closures.h"
#include "RHS.h"

using namespace amrex;

namespace PROB {

// problem parameters
struct ProbParm {
  Real rho_l = 0.072;  // density  [kg/m^3]
  Real u_l = 0.0;      // velocity [m/s]
  Real p_l = 7173.0;   // pressure [Pa]
  GpuArray<Real, NUM_SPECIES> Y_l = {
      0., 0.01277243, 0.,         0., 0., 0.10136214, 0.,
      0., 0.,         0.88586543, 0., 0., 0.};  // mass fractions [-] (molefrac 2:1:7)

  Real rho_r = 0.18075;
  Real u_r = -487.34;
  Real p_r = 35594.0;
  GpuArray<Real, NUM_SPECIES> Y_r = {
      0., 0.01277243, 0.,         0., 0., 0.10136214, 0.,
      0., 0.,         0.88586543, 0., 0., 0.};
};

inline Vector<std::string> cons_vars_names = indicies_t::get_cons_vars_names();
inline Vector<int> cons_vars_type = indicies_t::get_cons_vars_type();

typedef closures_dt<indicies_t, visc_suth_t, cond_suth_t,
                    multispecies_gas_t<indicies_t>> ProbClosures;
typedef rhs_dt<riemann_t<false, ProbClosures>, no_diffusive_t,
               reactor_t<ProbClosures>> ProbRHS;


void inline inputs() {}

// initial condition
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void prob_initdata(
    int i, int j, int k, Array4<Real> const &state,
    GeometryData const &geomdata, ProbClosures const &cls,
    ProbParm const &prob_parm) {
  const Real *prob_lo = geomdata.ProbLo();
  const Real *prob_hi = geomdata.ProbHi();
  const Real *dx = geomdata.CellSize();

  Real x = prob_lo[0] + (i + Real(0.5)) * dx[0];
  Real Pt, rhot, uxt;
  const Real* Yt;
  if (x < prob_hi[0] / 2) {
    Pt = prob_parm.p_l;
    rhot = prob_parm.rho_l;
    uxt = prob_parm.u_l;
    Yt = prob_parm.Y_l.data();
  } else {
    Pt = prob_parm.p_r;
    rhot = prob_parm.rho_r;
    uxt = prob_parm.u_r;
    Yt = prob_parm.Y_r.data();
  }
  Real et;
  cls.RYP2E(rhot, Yt, Pt, et);

  state(i, j, k, cls.UMX) = rhot * uxt;
  state(i, j, k, cls.UMY) = Real(0.0);
  state(i, j, k, cls.UMZ) = Real(0.0);
  state(i, j, k, cls.UET) = rhot * et + Real(0.5) * rhot * uxt * uxt;
  for (int n = 0; n < NUM_SPECIES; ++n) {
    state(i, j, k, cls.UFS + n) = rhot * Yt[n];
  }
}

// boundary conditions
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void bcnormal(
    const Real x[AMREX_SPACEDIM], Real dratio, const Real s_int[5],
    const Real s_refl[ProbClosures::NCONS], Real s_ext[5], const int idir,
    const int sgn, const Real time, GeometryData const & /*geomdata*/,
    ProbClosures const &closures, ProbParm const &prob_parm) {
  if (idir == 1) {  // ylo or yhi

    Abort("bcnormal not coded");   
  }
}

// source term
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void user_source(
    int i, int j, int k, const auto &state, const auto &rhs,
    const ProbParm &lprobparm, ProbClosures const &closures, auto const dx) {}
////////////////////////////////////////////////////////////////////////////////

///////////////////////////////AMR//////////////////////////////////////////////
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void user_tagging(
    int i, int j, int k, int nt_level, auto &tagfab, const auto &sdatafab,
    const auto &geomdata, const ProbParm &prob_parm, int level) {
}
////////////////////////////////////////////////////////////////////////////////

}  // namespace PROB
#endif