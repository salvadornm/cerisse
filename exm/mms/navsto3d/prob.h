#ifndef CNS_PROB_H_
#define CNS_PROB_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_ParmParse.H>
#include <Closures.h>
#include <RHS.h>
#include <Constants.h>
#include "mms.h"

using namespace amrex;

namespace PROB {

  static constexpr Real Pr   = 0.72;
  static constexpr Real gam  = 1.4;
  static constexpr Real mmw  = 29.0/1000.0;
  static constexpr Real Rgas = gas_constant/mmw;
  static constexpr Real Cp = Rgas*gam/(gam - 1.0);
  
  // problem parameters
  struct ProbParm {
    static constexpr Real u0= 1.0;
    static constexpr Real p0= 1.0;
    static constexpr Real rho0= 1.0;
  };

  struct viscparm_t {

    static constexpr int order   = 6;
    static constexpr bool use_LES = false;
    
    static constexpr Real viscosity     = 176.32;            // from python file (to get Re=1)
    static constexpr Real conductivity  = 245738.561825418;  // from python file (and Pr=0.72)
  };

  struct gasparm_t {
    static constexpr Real gamma            = gam;
    static constexpr Real molecular_weight = mmw;
  };

  // numerical method parameters
  struct methodparm_t {

    public:

    static constexpr bool dissipation = false;         // no dissipation
    static constexpr int  order = 2;                  // order numerical scheme
    static constexpr Real C2skew=0.1,C4skew=0.0016;   // Skew symmetric default

  };

typedef closures_dt<indicies_t, transport_const_t<viscparm_t>,
                    perfect_gas_t<gasparm_t, indicies_t> > ProbClosures;

template <typename cls_t > class user_source_t;

// HLLC-Riemann MUSCL
//typedef rhs_dt<riemann_t<false, ProbClosures>, viscous_t<viscparm_t,ProbClosures>, user_source_t <ProbClosures> >  ProbRHS;
// Skew
//typedef rhs_dt<skew_t<methodparm_t, ProbClosures>, viscous_t<viscparm_t,ProbClosures>, user_source_t <ProbClosures> > ProbRHS;
// Rusanov
//typedef rhs_dt<rusanov_t<ProbClosures>, viscous_t<viscparm_t,ProbClosures>, user_source_t <ProbClosures> >  ProbRHS;
// WENO & TENO   WenoZ5/Teno5/Teno6
//typedef rhs_dt<weno_t<ReconScheme::WenoZ5, ProbClosures>, viscous_t<viscparm_t,ProbClosures>, user_source_t <ProbClosures> > ProbRHS;
// KEEP 2/4/6
//typedef rhs_dt<keep_euler_t<false,false,4, ProbClosures>, viscous_t<viscparm_t,ProbClosures>, user_source_t <ProbClosures> > ProbRHS;
// CD 2/4/6
typedef rhs_dt<centraldif_t<false,false,6, ProbClosures>, viscous_t<viscparm_t,ProbClosures>, user_source_t <ProbClosures> > ProbRHS;

void inline inputs() {
  amrex::Print() << " MMS Euler " << std::endl;
}

// initial condition
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
prob_initdata(int i, int j, int k, Array4<Real> const &state,
              GeometryData const &geomdata, ProbClosures const &cls,
              ProbParm const &prob_parm) {
  
                const Real *prob_lo = geomdata.ProbLo();
  const Real *dx = geomdata.CellSize();
  Real x = prob_lo[0] + (i + Real(0.5)) * dx[0];
  Real y = prob_lo[1] + (j + Real(0.5)) * dx[1];
  Real z = prob_lo[2] + (k + Real(0.5)) * dx[2];
  
  const Real pi = 3.14159265358979323846; // in C++20 std::numbers::pi

  // initial conditions from mms
  Real rhot,ut,vt,wt,Pt;
  mms_exact(x, y, z, rhot, ut, vt, wt, Pt);

  
  state(i, j, k, cls.URHO) = rhot;
  state(i, j, k, cls.UMX)  = rhot * ut;
  state(i, j, k, cls.UMY)  = rhot * vt;
  state(i, j, k, cls.UMZ)  = rhot * wt;
  Real et = Pt / (cls.gamma - Real(1.0));
  state(i, j, k, cls.UET) = et + Real(0.5) * rhot * (ut * ut + vt*vt + wt*wt);
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
bcnormal(const Real x[AMREX_SPACEDIM], Real dratio, const Real s_int[5],
         const Real s_refl[ProbClosures::NCONS], Real s_ext[5], const int idir,
         const int sgn, const Real time, GeometryData const & /*geomdata*/,
         ProbClosures const &closures, ProbParm const &prob_parm) {
  if (idir == 1) { // ylo or yhi

    Abort("bcnormal not coded");
  }
}

// source term
///////////////////////////////SOURCE TERM /////////////////////////////////////
template <typename cls_t>
class user_source_t {
  public:
  void inline src(const Geometry& geomdata, const amrex::MFIter &mfi,
                  const amrex::Array4<const amrex::Real> &prims,
                  const amrex::Array4<amrex::Real> &rhs, const cls_t *cls_d,
                  amrex::Real dt){

    const Box bx = mfi.tilebox();

    // get geomdata!
    const Real *prob_lo = geomdata.ProbLo();
    const Real *dx = geomdata.CellSize();

    amrex::ParallelFor(bx,
      [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {

        Real x = prob_lo[0] + (i + Real(0.5)) * dx[0];
        Real y = prob_lo[1] + (j + Real(0.5)) * dx[1];
        Real z = prob_lo[2] + (k + Real(0.5)) * dx[2];

        // source terms from mms
        Real Srho,Srhou,Srhov,Srhow,Srhoe;
        mms_source(x, y, z, Srho,Srhou,Srhov,Srhow,Srhoe);
      
        //  MMS Source           
        rhs(i,j,k,cls_t::URHO) += Srho;
        rhs(i,j,k,cls_t::UMX)  += Srhou;
        rhs(i,j,k,cls_t::UMY)  += Srhov;
        rhs(i,j,k,cls_t::UMZ)  += Srhow;
        rhs(i,j,k,cls_t::UET)  += Srhoe;
       });
  };
};

///////////////////////////////AMR//////////////////////////////////////////////
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
user_tagging(int i, int j, int k, int nt_level, auto &tagfab,
             const auto &sdatafab, const auto &geomdata,
             const ProbParm &prob_parm, int level) {
}
////////////////////////////////////////////////////////////////////////////////

} // namespace PROB
#endif
