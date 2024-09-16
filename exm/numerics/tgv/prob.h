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
  Real Re = 1600.0;    
  Real Ma =  0.1;
  Real Pr = 0.71;

  Real omega_x = 1.0;
  Real omega_y = 1.0;
  Real omega_z = 1.0;

  Real L  = 1.0;
  Real v0 = 1.0;

};

inline Vector<std::string> cons_vars_names={"Xmom","Ymom","Zmom","Energy","Density"};
inline Vector<int> cons_vars_type={1,2,3,0,0};

typedef closures_dt<indicies_t, visc_suth_t, cond_suth_t,
                    calorifically_perfect_gas_t<indicies_t>>
    ProbClosures;

template <typename cls_t > class user_source_t;


typedef rhs_dt<riemann_t<false, ProbClosures>, no_diffusive_t, no_source_t > ProbRHS;
//typedef rhs_dt<skew_t<true,false, 4, ProbClosures>, no_diffusive_t, no_source_t > ProbRHS;


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

  const Real x = prob_lo[0] + (i + Real(0.5)) * dx[0];
  const Real y = prob_lo[1] + (j + Real(0.5)) * dx[1];
  const Real z = prob_lo[2] + (k + Real(0.5)) * dx[2];


  Real u[3], T, P,rhot,eint;

  const Real gam_M2 = prob_param.Ma*prob_param.Ma*cls.gamma;
  const Real p0 = 1.0/gam_M2;

  // TGV 
  // u
  u[0] = prob_parm.v0 * sin(prob_parm.omega_x * x / prob_parm.L) *
         cos(prob_parm.omega_y * y / prob_parm.L) *
         cos(prob_parm.omega_z * z / prob_parm.L);
  u[1] = -prob_parm.v0 * cos(prob_parm.omega_x * x / prob_parm.L) *
         sin(prob_parm.omega_y * y / prob_parm.L) *
         cos(prob_parm.omega_z * z / prob_parm.L);
  u[2] = Real(0.0);

  // P
  p =   p0 + prob_parm.rho0 * prob_parm.v0 * prob_parm.v0 / 16.0 *
                     (cos(2.0 * prob_parm.omega_x * x / prob_parm.L) +
                      cos(2.0 * prob_parm.omega_y * y / prob_parm.L)) *
                     (cos(2.0 * prob_parm.omega_z * z / prob_parm.L) + 2.0);
  // T 
  T = 1;

  // eint and rho
  rhot = gam_M2*p*T;

  // final state
  state(i, j, k, cls.URHO) = rhot;
  state(i, j, k, cls.UMX)  = rhot * u[0];
  state(i, j, k, cls.UMY)  = rhot * u[1];
  state(i, j, k, cls.UMZ)  = rhot * u[2];
  kin = Real(0.5)*(u[0] * u[0] + u[1] * u[1] + u[2] * u[2]);
  state(i, j, k, cls.UET)  = rhot*(eint + kin); 

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
