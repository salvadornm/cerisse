#ifndef CNS_PROB_H_
#define CNS_PROB_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_ParmParse.H>
#include <Closures.h>
#include <RHS.h>
#include <bc_types.h>

#include <nscbc.h>

// 2D Convective Vortex 


using namespace amrex;

namespace PROB {

// 2-D Lodato-style outward pressure-pulse test.
// Quiet ideal-gas field, centered Gaussian pressure pulse, all faces outflow/NSCBC.
struct ProbParm {
  Real gamma = 1.4;
  Real p0    = 101325.0;
  Real T0    = 300.0;
  Real Rair  = 287.0;
  Real rho0  = p0 / (Rair * T0);
  Real c0    = std::sqrt(gamma * p0 / rho0);

  Real Lx    = 0.16;
  Real Ly    = 0.16;
  Real delta = 0.01;   // 1% pressure pulse, as in the 2-D version in Rathore thesis (bigger than Lodato)
  Real Rp    = 0.02;  // Rp  is 12/5% domain pulse radius for L=0.16 m (bigger than Lodato)

  Real Y0[NUM_SPECIES] = {1.0};
};


// numerical method parameters
struct methodparm_t {

  public:

  static constexpr bool dissipation = true;         // no dissipation
  static constexpr int  order = 6;                  // order numerical scheme   
  static constexpr Real C2skew=0.1,C4skew=0.016;   // Skew symmetric default

};

typedef closures_dt<indicies_t, visc_suth_t, cond_suth_t,
                    calorifically_perfect_gas_t<indicies_t>  > ProbClosures;

//template <typename cls_t > class user_source_t;


//typedef rhs_dt<riemann_t<false, ProbClosures>, no_diffusive_t, no_source_t > ProbRHS;
//typedef rhs_dt<skew_t<methodparm_t, ProbClosures>, no_diffusive_t, no_source_t > ProbRHS;
typedef rhs_dt<weno_t<ReconScheme::Teno5, ProbClosures>, no_diffusive_t,  no_source_t>  ProbRHS;


// boundary conditions
typedef manual_bc_t<ProbClosures> GlobalBC;


void inline inputs() {
  ProbParm data;
  amrex::Print() << "************** 2-D pressure-pulse NSCBC test **************\n";
  amrex::Print() << "p0=" << data.p0 << " Pa, T0=" << data.T0
                 << " K, c0=" << data.c0 << " m/s\n";
  amrex::Print() << "Lx=" << data.Lx << ", Ly=" << data.Ly
                 << ", delta=" << data.delta << ", Rp=" << data.Rp << "\n";
  amrex::Print() << "All physical faces should be cns.*_bc=2 and cns.nscbc_*=2.\n";
  amrex::Print() << "***********************************************************\n";
}


// initial condition
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
prob_initdata(int i, int j, int k, Array4<Real> const &state,
              GeometryData const &geomdata, ProbClosures const &cls,
              ProbParm const &prob_parm) {
  const Real *prob_lo = geomdata.ProbLo();
  const Real *prob_hi = geomdata.ProbHi();
  const Real *dx = geomdata.CellSize();

  const Real x = prob_lo[0] + (Real(i) + Real(0.5)) * dx[0];
  const Real y = prob_lo[1] + (Real(j) + Real(0.5)) * dx[1];

  const Real xc = Real(0.0); //Real(0.5) * (geomdata.ProbLo(0) + geomdata.ProbHi(0));
  const Real yc = Real(0.0); //Real(0.5) * (geomdata.ProbLo(1) + geomdata.ProbHi(1));

  const Real r2 = (x - xc) * (x - xc) + (y - yc) * (y - yc);
  
  const Real P  = prob_parm.p0 *
                  (Real(1.0) + prob_parm.delta *
                   std::exp(-r2 / (Real(2.0) * prob_parm.Rp * prob_parm.Rp)));

  // Isothermal initialization, exactly as Lodato's pressure-pulse test: rho=p/(R T0).
  const Real rho = P / (prob_parm.Rair * prob_parm.T0);
  const Real u = Real(0.0), v = Real(0.0), w = Real(0.0);
  const Real rhoeint = P / (prob_parm.gamma - Real(1.0));

#if NUM_SPECIES > 1
  for (int n = 0; n < NUM_SPECIES; ++n) {
    state(i,j,k,cls.UFS+n) = rho * prob_parm.Y0[n];
  }
#else
  state(i,j,k,cls.URHO) = rho;
#endif
  state(i,j,k,cls.UMX) = rho * u;
  state(i,j,k,cls.UMY) = rho * v;
  state(i,j,k,cls.UMZ) = rho * w;
  state(i,j,k,cls.UET) = rhoeint + Real(0.5) * rho * (u*u + v*v + w*w);
}

/////////////////////////////// BC /////////////////////////////////////////////
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
bcnormal(const Real x[AMREX_SPACEDIM], Real dratio, const Real s_int[ProbClosures::NCONS],
         const Real s_refl[ProbClosures::NCONS], Real s_ext[ProbClosures::NCONS], const int idir,
         const int sgn, const Real time, GeometryData const & /*geomdata*/,
         ProbClosures const &closures, ProbParm const &prob_parm) {

  //        
  for (int n = 0; n < ProbClosures::NCONS; ++n) {
    s_ext[n] = s_int[n];
  }

  return;
  //

  const int face = (idir+1)*sgn;

  switch(face)
  {
    case  2:  // SOUTH
      break;
    case  1:  // WEST x=0
      //GlobalBC::bc_inlet_fixmassflow(1.0,0.0,0.0,&closures,
      //      prob_parm.Q,prob_parm.T0,prob_parm.Y0, s_int, s_ext);
      
      break;
    case -1:  // EAST x= Lx
     // GlobalBC::bc_fixP(-1.0,0.0,0.0,&closures,prob_parm.p0, s_int, s_ext);
     //GlobalBC::bc_subsonic_outflow_fixP(-1.0,0.0,0.0,&closures,prob_parm.p0, s_int, s_ext);
     
      break;
    case -2:   // NORTH
      break;
    default:

      break;
  }
}
///////////////////////////////AMR//////////////////////////////////////////////
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
user_tagging(int i, int j, int k, int nt_level, auto &tagfab,
             const auto &sdatafab, const auto &geomdata,
             const ProbParm &prob_parm, int level) {

}
///////////////////////////////////////////////////////////////////////////////

} // namespace PROB
#endif
