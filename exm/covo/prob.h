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
  Real mach = 0.05;
  Real beta = 0.02; // vortex strength
  Real p0 = 1e6;    // [erg cm^-3]
  Real T0 = 300.0;  // [K]
  Real rho0;
  Real v0 = 1735.95;
  // v0 = 0.05 * np.sqrt(1.4 * 287 * 300) * 100 # convection velocity 1735.95
};

inline Vector<std::string> cons_vars_names={"Xmom","Ymom","Zmom","Energy","Density"};
inline Vector<int> cons_vars_type={1,2,3,0,0};

typedef closures_dt<indicies_t, visc_suth_t, cond_suth_t,
                    calorifically_perfect_gas_t<indicies_t>>
    ProbClosures;

template <typename cls_t > class user_source_t;



typedef rhs_dt<riemann_t<false, ProbClosures>, no_diffusive_t, no_source_t > ProbRHS;


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
  
  // Vortex functions 
  const Real xc = 5.0;
  const Real yc = 5.0; // vortex initial pos
  const Real R = 0.5;  // vortex radius
  const Real rsq = (x - xc) * (x - xc) + (y - yc) * (y - yc);
  amrex::Real  cp  = 1000; // temp snm

  amrex::Real u[3], T, Pt,rhot;

  // auxiliar parameters
  const Real vbeta = prob_parm.v0 * prob_parm.beta;
  const Real expd = exp(-0.5 * rsq / R / R);
  

  // velocity
  u[0] = prob_parm.v0  - vbeta * (y - yc) / R * expd;
  u[1] = vbeta * (x - xc) / R * expd;
  u[2] = 0.0;
  // Temperature
  T = prob_parm.T0 - Real(0.5) * (vbeta*vbeta) / cp * exp(-rsq / R / R);
  // Pressure

  // Density  and Internal Energy
    
  Real eint = Pt / (cls.gamma - Real(1.0));

  // final state
  state(i, j, k, cls.URHO) = rhot;
  state(i, j, k, cls.UMX)  = rhot * u[0];
  state(i, j, k, cls.UMY)  = rhot * u[1];
  state(i, j, k, cls.UMZ)  = Real(0.0);  
  state(i, j, k, cls.UET) = rhot*(eint + Real(0.5) * (u[0] * u[0] + u[1] * u[1])); 

}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
bcnormal(const Real x[AMREX_SPACEDIM], Real dratio, const Real s_int[5],
         const Real s_refl[ProbClosures::NCONS], Real s_ext[5], const int idir,
         const int sgn, const Real time, GeometryData const & /*geomdata*/,
         ProbClosures const &closures, ProbParm const &prob_parm) {
}

///////////////////////////////AMR//////////////////////////////////////////////
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
user_tagging(int i, int j, int k, int nt_level, auto &tagfab,
             const auto &sdatafab, const auto &geomdata,
             const ProbParm &prob_parm, int level) {


  Real rhot = sdatafab(i,j,k,ProbClosures::URHO);

  Real dengrad_threshold = 0.5;
  Real drhox = Math::abs(sdatafab(i+1,j,k,ProbClosures::URHO) -
   sdatafab(i,j,k,ProbClosures::URHO))/rhot;

  Real drhoy = Math::abs(sdatafab(i,j+1,k,ProbClosures::URHO) -
   sdatafab(i,j-1,k,ProbClosures::URHO))/rhot;

  if (nt_level > 0) 
  {
    // tag cells based on density values
    switch (level)
    {
      case 0:
        tagfab(i,j,k) = (rhot > 1.1 && rhot < 1.9);
        break;
      case 1:
        tagfab(i,j,k) = (rhot > 1.2 && rhot < 1.8); 
        break;
      default:
        tagfab(i,j,k) = (rhot > 1.3 && rhot < 1.7); 
        break;
    }
  }
} 
////////////////////////////////////////////////////////////////////////////////

} // namespace PROB
#endif
