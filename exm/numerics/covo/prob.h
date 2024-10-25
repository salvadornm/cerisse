#ifndef CNS_PROB_H_
#define CNS_PROB_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_ParmParse.H>
#include <Closures.h>
#include <RHS.h>

// 2D Convective Vortex 


using namespace amrex;

namespace PROB {

// problem parameters  (1:bottom   2:top)
struct ProbParm {
  Real p0 = 101300.0;    
  Real rho0 =  1.17170407;
  Real v0 = 35.0;
  Real beta = 0.04; // parameter
  Real L = 0.3112;  // size of domain
};

inline Vector<std::string> cons_vars_names={"Xmom","Ymom","Zmom","Energy","Density"};
inline Vector<int> cons_vars_type={1,2,3,0,0};

typedef closures_dt<indicies_t, visc_suth_t, cond_suth_t,
                    calorifically_perfect_gas_t<indicies_t>>
    ProbClosures;

template <typename cls_t > class user_source_t;


//typedef rhs_dt<riemann_t<false, ProbClosures>, no_diffusive_t, no_source_t > ProbRHS;
typedef rhs_dt<skew_t<false,false, 4, ProbClosures>, no_diffusive_t, no_source_t > ProbRHS;
// typedef rhs_dt<weno_t<ReconScheme::Teno5, ProbClosures>, no_diffusive_t, no_source_t > ProbRHS;


void inline inputs() {
  ParmParse pp;

  pp.add("cns.order_rk", 3);   // -2, 1, 2 or 3"
  pp.add("cns.stages_rk", 3);  // 1, 2 or 3
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
  const Real xc = 0.5*prob_parm.L; const Real yc = 0.5*prob_parm.L;
  const Real R  = 0.05*prob_parm.L; // vortex radius  L/20
  const Real rsq = (x - xc) * (x - xc) + (y - yc) * (y - yc);

  amrex::Real u[3], T, P,rhot;

  // auxiliar parameters
  const Real expd = exp(- rsq / R / R);
  const Real Gam  = prob_parm.beta*prob_parm.v0*R*sqrt(exp(Real(1.0)));
  const Real Pfluc  = 2.0 *prob_parm.rho0*(Gam/R)*(Gam/R);

  
  // psi = Gamm*expd    stream function
  // vx = dpsi/dy  and vy =  -dpsi/dx

  // velocity
  if (rsq < 4*R) 
  {
    u[0] = prob_parm.v0  - Real(2.0)*Gam/(R*R) * expd * (y - yc);
    u[1] = Real(2.0)*Gam/(R*R) * expd * (x - xc);
    u[2] = Real(0.0);
  }
  else
  {
    u[0]=0.0; u[1]=0.0;u[2]=0.0;	  
  }
  // Pressure
  P = prob_parm.p0 - Pfluc*expd;

  // Density  and Internal Energy

  rhot = prob_parm.rho0;  
  Real eint = P / (cls.gamma - Real(1.0));

  // final state
  state(i, j, k, cls.URHO) = rhot;
  state(i, j, k, cls.UMX)  = rhot * u[0];
  state(i, j, k, cls.UMY)  = rhot * u[1];
  state(i, j, k, cls.UMZ)  = Real(0.0);  
  state(i, j, k, cls.UET)  = eint + Real(0.5) * rhot * (u[0] * u[0] + u[1] * u[1]); 

}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
bcnormal(const Real x[AMREX_SPACEDIM], Real dratio, const Real s_int[5],
         const Real s_refl[ProbClosures::NCONS], Real s_ext[5], const int idir,
         const int sgn, const Real time, GeometryData const & /*geomdata*/,
         ProbClosures const &closures, ProbParm const &prob_parm) {
}
//
///////////////////////////////AMR//////////////////////////////////////////////
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
user_tagging(int i, int j, int k, int nt_level, auto &tagfab,
             const auto &sdatafab, const auto &geomdata,
             const ProbParm &prob_parm, int level) {

}
///////////////////////////////////////////////////////////////////////////////

} // namespace PROB
#endif
