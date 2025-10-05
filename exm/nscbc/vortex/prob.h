#ifndef CNS_PROB_H_
#define CNS_PROB_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_ParmParse.H>
#include <Closures.h>
#include <RHS.h>
#include <bc_types.h>

// 2D Convective Vortex 


using namespace amrex;

namespace PROB {

// problem parameters  (1:bottom   2:top)
struct ProbParm {
  Real gamma = 1.4;        // ratio of specific heats
  Real p0 = 101325.0;      
  Real T0 = 300.0;  
  Real Rair = 287.0;                 // specific gas constant
  Real rho0 = p0 / (Rair * T0);      // Ideal gas law
  Real c0 = sqrt(gamma * p0 / rho0); // speed of soundR
  Real Ma = 0.575;
  Real u0 = Ma * c0;       // reference velocity
  Real Lx = 1.0;           // domain length in x
  Real Ly = 1.0;           // domain length in y  
  Real beta = 0.04;        // vortex strength
  Real Y0[NUM_SPECIES] = {0.0};
  Real Q =  rho0*u0;      // incoming flow rate (per area)
};


// numerical method parameters
struct methodparm_t {

  public:

  static constexpr bool dissipation = false;         // no dissipation
  static constexpr int  order = 4;                  // order numerical scheme   
  static constexpr Real C2skew=0.1,C4skew=0.16;   // Skew symmetric default

};


typedef closures_dt<indicies_t, visc_suth_t, cond_suth_t,
                    calorifically_perfect_gas_t<indicies_t> > ProbClosures;

//template <typename cls_t > class user_source_t;


typedef rhs_dt<riemann_t<false, ProbClosures>, no_diffusive_t, no_source_t > ProbRHS;
//typedef rhs_dt<skew_t<methodparm_t, ProbClosures>, no_diffusive_t, no_source_t > ProbRHS;

// boundary conditions
typedef manual_bc_t<ProbClosures> GlobalBC;


void inline inputs() {
  ProbParm data;
  amrex::Print() << "**************                          ************ " << std::endl;
  amrex::Print() << " Two-dimensional vortex convection " << std::endl;
  amrex::Print() << "**************                          ************ " << std::endl;
  Real Ma = data.u0 / data.c0;
  amrex::Print() << " Mach   =  " << Ma << std::endl;
  amrex::Print() << "**************                          ************ " << std::endl;
}

// initial condition
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
prob_initdata(int i, int j, int k, Array4<Real> const &state,
              GeometryData const &geomdata, ProbClosures const &cls,
              ProbParm const &prob_parm) {
  const Real *prob_lo = geomdata.ProbLo();
  const Real *prob_hi = geomdata.ProbHi();
  const Real *dx = geomdata.CellSize();

  Real x = prob_lo[0] + (i + Real(0.5)) * dx[0];
  Real y = prob_lo[1] + (j + Real(0.5)) * dx[1];
  
  // Vortex position (xc,yc) middle of domain
  const Real xc = 0.5*prob_parm.Lx; const Real yc = 0.5*prob_parm.Ly;
  const Real Rvtx  = 0.05*prob_parm.Lx; // vortex radius  L/20
  const Real rsq = (x - xc) * (x - xc) + (y - yc) * (y - yc);

  amrex::Real u[3]={0.0}, T, P,rhot;

  // auxiliar parameters
  const Real expd   = exp(- rsq / Rvtx / Rvtx);
  const Real Gam    = prob_parm.beta*prob_parm.u0*Rvtx*sqrt(exp(Real(1.0)));
  const Real Pfluc  = 2.0 *prob_parm.rho0*(Gam/Rvtx)*(Gam/Rvtx);

  // main Velocity
  u[0] = prob_parm.u0;
  // main Pressure
  P = prob_parm.p0;
  
  //------
  // psi = Gamm*expd    stream function
  // vx = dpsi/dy  and vy =  -dpsi/dx

  // vortex only defined for r < 2*Rvtx 
  if (rsq < 4.0*Rvtx) 
  {
    u[0] += - Real(2.0)*Gam/(Rvtx*Rvtx) * expd * (y - yc);
    u[1] +=   Real(2.0)*Gam/(Rvtx*Rvtx) * expd * (x - xc);
    P    +=   - Pfluc*expd;
     //P += 0.5*Pfluc*exp(-2.0*rsq / Rvtx / Rvtx);
  }  
  
  //----------------

  // Density  
  //cls.PYT2R(P, prob_parm.Y0, prob_parm.T0, rhot);
  rhot = prob_parm.rho0;
  // Internal energy  (rhoe = P/(gamma-1)
  Real eint = P / (cls.gamma - Real(1.0));

  // final state
  state(i, j, k, cls.URHO) = rhot;
  state(i, j, k, cls.UMX)  = rhot * u[0];
  state(i, j, k, cls.UMY)  = rhot * u[1];
  state(i, j, k, cls.UMZ)  = Real(0.0);  
  state(i, j, k, cls.UET)  = eint + Real(0.5) * rhot * (u[0] * u[0] + u[1] * u[1]); 

}

/////////////////////////////// BC /////////////////////////////////////////////
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
bcnormal(const Real x[AMREX_SPACEDIM], Real dratio, const Real s_int[ProbClosures::NCONS],
         const Real s_refl[ProbClosures::NCONS], Real s_ext[ProbClosures::NCONS], const int idir,
         const int sgn, const Real time, GeometryData const & /*geomdata*/,
         ProbClosures const &closures, ProbParm const &prob_parm) {

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
      //GlobalBC::bc_fixP(-1.0,0.0,0.0,&closures,prob_parm.p0, s_int, s_ext);
      GlobalBC::bc_subsonic_outflow_fixP(-1.0,0.0,0.0,&closures,prob_parm.p0, s_int, s_ext);
     
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
