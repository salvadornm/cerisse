#ifndef CNS_PROB_H_
#define CNS_PROB_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_ParmParse.H>
#include <PelePhysics.H>

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
      0.01736973, 0.13784636, 0., 0.,        0.,
      0.,         0.,         0., 0.84478391};  // mass fractions [-]

  Real rho_r = 0.18075;
  Real u_r = -487.34;
  Real p_r = 35594.0;
  GpuArray<Real, NUM_SPECIES> Y_r = {0.01736973, 0.13784636, 0., 0.,        0.,
                                     0.,         0.,         0., 0.84478391};
};

// these two are better to keep in setup/index?
inline Vector<std::string> cons_vars_names = indicies_t::get_cons_vars_names();
inline Vector<int> cons_vars_type = indicies_t::get_cons_vars_type();

typedef closures_dt<indicies_t, visc_suth_t, cond_suth_t,
                    multispecies_gas_t<indicies_t>>
    ProbClosures;
typedef rhs_dt<riemann_t<false, ProbClosures>, no_diffusive_t, no_source_t>
    ProbRHS;

void inline inputs() {
  ParmParse pp;

  pp.add("cns.order_rk", 3);   // -2, 1, 2 or 3"
  pp.add("cns.stages_rk", 3);  // 1, 2 or 3
}

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
  state(i, j, k, cls.UMX) = rhot * uxt;
  state(i, j, k, cls.UMY) = Real(0.0);
  state(i, j, k, cls.UMZ) = Real(0.0);
  Real et;
  cls.RYP2E(rhot, Yt, Pt, et);
  state(i, j, k, cls.UET) = rhot * et + Real(0.5) * rhot * uxt * uxt;
  for (int n = 0; n < NUM_SPECIES; ++n) {
    state(i, j, k, cls.UFS + n) = rhot * Yt[n];
  }
}

// boundary conditions
/**
 * \brief Fill external boundary conditions for ghost cells.
 *
 * @param x         ghost cell cooridinates.
 * @param dr        wall-ghost/wall-first internal distance ratio
 * @param s_int     flow state inside of the domain.
 * @param s_ext     flow state to be filled.
 * @param idir      direction (0: x, 1: y, 2: z).
 * @param sgn       high or low boundary (1: low, -1: high).
 * @param time      time.
 * @param geomdata  domain geometry data.
 * @param prob_parm ProbParm data as defined in prob_parm.H and initialised in
 * amrex_probinit.
 * @sa CnsFillExtDir
 * @sa CnsFillExtDir::operator()
 */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void bcnormal(
    const Real x[AMREX_SPACEDIM], Real dratio, const Real s_int[5],
    const Real s_refl[ProbClosures::NCONS], Real s_ext[5], const int idir,
    const int sgn, const Real time, GeometryData const & /*geomdata*/,
    ProbClosures const &closures, ProbParm const &prob_parm) {
  if (idir == 1) {  // ylo or yhi

    Abort("bcnormal not coded");

    // Real q_ext[NPRIM] = {0.0};
    // no-slip
    // q_ext[QU]    = -s_int[UMX]/s_int[URHO];
    // q_ext[QV]    = -s_int[UMY]/s_int[URHO];
    // q_ext[QW]    = -s_int[UMZ]/s_int[URHO];

    // // dp/dn = 0
    // Real eint_int = (s_int[UET] - 0.5*(s_int[UMX]*s_int[UMX] +
    // s_int[UMY]*s_int[UMY] + s_int[UMZ]*s_int[UMZ])/s_int[URHO])/s_int[URHO];
    // Real p_int = (parm.eos_gamma - 1.0)*s_int[URHO]*eint_int;
    // q_ext[QPRES] = p_int;
    // // T=Twall
    // Real T_int = p_int/(parm.Rspec*s_int[URHO]);
    // q_ext[QT]    = max(prob_parm.Tw  +  dratio*(prob_parm.Tw - T_int),50.0);
    // // rho = eos(P,T)
    // q_ext[QRHO]  = q_ext[QPRES]/(parm.Rspec*q_ext[QT]);

    // // convert prims to cons
    // s_ext[URHO] = q_ext[QRHO];
    // s_ext[UMX] = q_ext[QRHO]*q_ext[QU];
    // s_ext[UMY] = q_ext[QRHO]*q_ext[QV];
    // s_ext[UMZ] = q_ext[QRHO]*q_ext[QW];
    // Real ekin_ext = 0.5*(q_ext[QU]*q_ext[QU] + q_ext[QV]*q_ext[QV] +
    // q_ext[QW]*q_ext[QW]); Real eint_ext =
    // q_ext[QPRES]/(q_ext[QRHO]*(parm.eos_gamma - 1.0)); s_ext[UET] =
    // q_ext[QRHO]*(eint_ext + ekin_ext);
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
  // Real dengrad_threshold = 0.5;
  // Real drhox = Math::abs(sdatafab(i+1,j,k,URHO) -
  // sdatafab(i-1,j,k,URHO))/sdatafab(i,j,k,URHO); if (drhox >
  // dengrad_threshold) {
  //   tagfab(i,j,k) = true;
  //   tagfab(i+1,j,k) = true;
  //   tagfab(i+2,j,k) = true;
  //   tagfab(i+3,j,k) = true;
  // }
}
////////////////////////////////////////////////////////////////////////////////

}  // namespace PROB
#endif