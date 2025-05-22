#ifndef CNS_PROB_H_
#define CNS_PROB_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_ParmParse.H>
#include <Closures.h>
#include <RHS.h>
#include "Utilities.h"

using namespace amrex;

namespace PROB {

// problem parameters
struct ProbParm {
  Real u0= 1.0;
  Real p0= 100000.0;
  Real rho0= 1.0;
};

// numerical method parameters
struct methodparm_t {

  public:

  static constexpr bool dissipation = false;         // no dissipation
  static constexpr int  order = 2;                  // order numerical scheme
  static constexpr Real C2skew=0.1,C4skew=0.0016;   // Skew symmetric default

};

typedef closures_dt<indicies_t, visc_suth_t, cond_suth_t,
                    calorifically_perfect_gas_t<indicies_t>> ProbClosures;

template <typename cls_t > class user_source_t;

// HLLC-Riemann MUSCL
//typedef rhs_dt<riemann_t<false, ProbClosures>, no_diffusive_t,  no_source_t >  ProbRHS;
// Skew
//typedef rhs_dt<skew_t<methodparm_t, ProbClosures>, no_diffusive_t,  no_source_t > ProbRHS;
// Rusanov
//typedef rhs_dt<rusanov_t<ProbClosures>, no_diffusive_t,  no_source_t >  ProbRHS;
// WENO & TENO   WenoZ5/Teno5/Teno6
//typedef rhs_dt<weno_t<ReconScheme::WenoZ5, ProbClosures>, no_diffusive_t,  no_source_t> > ProbRHS;
// KEEP 2/4/6
//typedef rhs_dt<keep_euler_t<false,false,4, ProbClosures>, no_diffusive_t,  no_source_t> ProbRHS;
// CD 2/4/6
typedef rhs_dt<centraldif_t<false,false,4, ProbClosures>, no_diffusive_t, no_source_t > ProbRHS;

void inline inputs() {
  amrex::Print() << " HIT test  " << std::endl;
}

//-------------------------------------------------------------------------------------------
// initial condition
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void prob_initdata(
  int i, int j, int k, Array4<Real> const &state,
  GeometryData const &geomdata, ProbClosures const &cls,
  ProbParm const &prob_parm, Utility* util = nullptr) {

  const Real *prob_lo = geomdata.ProbLo();
  const Real *dx = geomdata.CellSize();
  Real x = prob_lo[0] + (i + Real(0.5)) * dx[0];
  Real y = prob_lo[1] + (j + Real(0.5)) * dx[1];
  Real z = prob_lo[2] + (k + Real(0.5)) * dx[2];
  
  // initial conditions rho P
  Real rhot = prob_parm.rho0;
  Real Pt   = prob_parm.p0;

  // read velocity field from datafile
  Real ut,vt,wt;
  util->get_velocity(i,j,k,ut,vt,wt);
  //

  state(i, j, k, cls.URHO) = rhot;
  state(i, j, k, cls.UMX)  = rhot * ut;
  state(i, j, k, cls.UMY)  = rhot * vt;
  state(i, j, k, cls.UMZ)  = rhot * wt;
  Real et = Pt / (cls.gamma - Real(1.0));
  state(i, j, k, cls.UET) = et + Real(0.5) * rhot * (ut * ut + vt*vt + wt*wt);
}
//-------------------------------------------------------------------------------------------
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
bcnormal(const Real x[AMREX_SPACEDIM], Real dratio, const Real s_int[5],
         const Real s_refl[ProbClosures::NCONS], Real s_ext[5], const int idir,
         const int sgn, const Real time, GeometryData const & /*geomdata*/,
         ProbClosures const &closures, ProbParm const &prob_parm) {
  if (idir == 1) { // ylo or yhi

    Abort("bcnormal not coded");
  }
}

///////////////////////////////AMR//////////////////////////////////////////////
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
user_tagging(int i, int j, int k, int nt_level, auto &tagfab,
             const auto &sdatafab, const auto &geomdata,
             const ProbParm &prob_parm, int level) {
}
////////////////////////////////////////////////////////////////////////////////

} // namespace PROB
#endif
