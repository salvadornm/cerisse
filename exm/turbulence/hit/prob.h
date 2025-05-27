#ifndef CNS_PROB_H_
#define CNS_PROB_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_ParmParse.H>

#include <Closures.h>
#include <RHS.h>
#include <LES.h>
#include "Utilities.h"

using namespace amrex;

namespace PROB {

// problem parameters
struct ProbParm {  
  Real rho0= 1.0;
  Real gamma = 1.4;
  Real L0 = 6.28318;
  Real T0 = 1.0;  
  Real u0 = 1.0;
  Real Ma_rms = 0.2; 
  Real c0 = u0/Ma_rms;
  Real p0 = rho0*c0*c0/gamma;
  Real t0 = L0/u0;
};

// numerical method parameters
struct methodparm_t {

  public:

  static constexpr bool dissipation = true;         // no dissipation
  static constexpr int  order = 4;                  // order numerical scheme
  static constexpr Real C2skew=0.1,C4skew=0.016;   // Skew symmetric default

};
// LES parameters
struct LESparm_t {

  public:

  static constexpr int  order = 2;         // order numerical scheme (for gradient estimation)
  static constexpr Real Pr_o_Prsgs = 0.8;  // Pr/Prsgs
  static constexpr Real Scsgs = 0.7;       // sgs Sc
  static constexpr Real Cs = 0.0;          // Smagorinsky constant
  static constexpr Real CI = 0.08;         // Yoshizawa constant
  static constexpr bool fixDelta = false;  // Fix filter witdth
  static constexpr Real Delta = 0.02;      // Filter width (if above true)  L/20
};


struct viscparm_t {

  public :

  static constexpr int  order = 2;
  static constexpr bool use_LES = true;

  // zero-viscosity and conductivity  (infinite Reynolds)
  static constexpr Real viscosity    = 0.0;  // viscosity    (for constant value)  
  static constexpr Real conductivity = 0.0;  // conductivity (for constant value)

};

// changes indicies_t to indicesgen_t<4>, which allocates more ghost points

typedef closures_dt<indicies_t, transport_const_t<viscparm_t>,
                    calorifically_perfect_gas_t<indicies_t>, Smagorinsky_t<LESparm_t,indicies_t> > ProbClosures;

template <typename cls_t > class user_source_t;

// HLLC-Riemann MUSCL
//typedef rhs_dt<riemann_t<false, ProbClosures>, no_diffusive_t,  no_source_t >  ProbRHS;
// Skew
typedef rhs_dt<skew_t<methodparm_t, ProbClosures>, no_diffusive_t,  no_source_t > ProbRHS;
// Rusanov
//typedef rhs_dt<rusanov_t<ProbClosures>, no_diffusive_t,  no_source_t >  ProbRHS;
// WENO & TENO   WenoZ5/Teno5/Teno6
//typedef rhs_dt<weno_t<ReconScheme::WenoZ5, ProbClosures>, no_diffusive_t,  no_source_t> > ProbRHS;
// KEEP 2/4/6
//typedef rhs_dt<keep_euler_t<false,false,4, ProbClosures>, no_diffusive_t,  no_source_t> ProbRHS;
// CD 2/4/6
//typedef rhs_dt<centraldif_t<false,false,4, ProbClosures>, no_diffusive_t, no_source_t > ProbRHS;

void inline inputs() {
  ProbParm data;
  amrex::Print() << "**************  " << std::endl;
  amrex::Print() << " HIT test  " << std::endl;
  amrex::Print() << " Ma_rms =  " << data.u0/data.c0 << std::endl;
  amrex::Print() << " c0 [m/s]=  " << data.c0 << std::endl;
  amrex::Print() << " u0 (RMS)[m/s]=  " << data.u0 << std::endl;
  amrex::Print() << " P0 [s]  =  " << data.p0 << std::endl;  
  amrex::Print() << " t0 [s]  =  " << data.t0 << std::endl;
  
  amrex::Print() << "**************  " << std::endl;
}

//-------------------------------------------------------------------------------------------
// initial condition
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void prob_initdata(
  int i, int j, int k, Array4<Real> const &state,
  GeometryData const &geomdata, ProbClosures const &cls,
  ProbParm const &prob_parm, Utility* util = nullptr) {

  // const Real *prob_lo = geomdata.ProbLo();
  // const Real *dx = geomdata.CellSize();
  // Real x = prob_lo[0] + (i + Real(0.5)) * dx[0];
  // Real y = prob_lo[1] + (j + Real(0.5)) * dx[1];
  // Real z = prob_lo[2] + (k + Real(0.5)) * dx[2];
  
  // initial conditions rho P
  Real rhot = prob_parm.rho0;
  Real Pt   = prob_parm.p0;
  Real u0   = prob_parm.u0;

  // read velocity field from datafile
  Real ut,vt,wt;
  util->get_velocity(i,j,k,ut,vt,wt);
  // scale to rms
  ut *= u0; vt *= u0; wt *= u0;
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
         const int /*sgn*/, const Real time, GeometryData const & /*geomdata*/,
         ProbClosures const &closures, ProbParm const &/*prob_parm*/) {
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
