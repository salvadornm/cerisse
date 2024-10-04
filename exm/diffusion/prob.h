#ifndef CNS_PROB_H_
#define CNS_PROB_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_ParmParse.H>
#include <Closures.h>
#include <RHS.h>
#include <Constants.H>

// 2D Thermal diffusion
// Pure diffusin case taken from HAMISH
// https://www.ukctrf.com/index.php/benchmarking-of-the-new-software/

using namespace amrex;

namespace PROB {

// problem parameters 
struct ProbParm {
  Real p = pres_atm2si;
  Real u = Real(0.0);
  Real v = Real(0.0);  
  Real T0 = Real(300.0);  
  Real T1 = Real(100.0);  
  
  Real L = Real(0.016); // domain size
  Real x0 = Real(0.5)*L;
  Real y0 = Real(0.5)*L;
  Real delta = Real(L/40.0);
  
};

inline Vector<std::string> cons_vars_names={"Xmom","Ymom","Zmom","Energy","Density"};
inline Vector<int> cons_vars_type={1,2,3,0,0};

typedef closures_dt<indicies_t, visc_suth_t, cond_suth_t,
                    calorifically_perfect_gas_t<indicies_t>> ProbClosures;

template <typename cls_t > class user_source_t;


typedef rhs_dt<no_euler_t, no_diffusive_t, no_source_t > ProbRHS;
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

  Real x = prob_lo[0] + (i + Real(0.5)) * dx[0];
  Real y = prob_lo[1] + (j + Real(0.5)) * dx[1];
  
  // local vars
  Real rho,eint, T, P,u[2]; 
  
  // initial conditioms
  u[0] = prob_parm.u;u[1] = prob_parm.v;
  P    = prob_parm.p;

  Real rad2   = (x - prob_parm.x0)*(x - prob_parm.x0) + (y - prob_parm.y0)*(y - prob_parm.y0);
  Real delta2 = prob_parm.delta*prob_parm.delta;
  T =  prob_parm.T0 + prob_parm.T1*std::exp(-0.25_rt*rad2/delta2);

  // rhp and E from P and T
  rho = P/(cls.Rspec*T);
  eint =  P / (cls.gamma - Real(1.0));

  state(i, j, k, cls.URHO) = rho;
  state(i, j, k, cls.UMX)  = rho * u[0];
  state(i, j, k, cls.UMY)  = rho * u[1];
  state(i, j, k, cls.UMZ)  = Real(0.0);  
  state(i, j, k, cls.UET)  = eint + Real(0.5) * rho * (u[0] * u[0] + u[1] * u[1]);    
}

/////////////////////////////// BC /////////////////////////////////////////////
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
bcnormal(const Real x[AMREX_SPACEDIM], Real dratio, const Real s_int[ProbClosures::NCONS],
         const Real s_refl[ProbClosures::NCONS], Real s_ext[ProbClosures::NCONS], const int idir,
         const int sgn, const Real time, GeometryData const & /*geomdata*/,
         ProbClosures const &closures, ProbParm const &prob_parm) {

  // const int URHO = ProbClosures::URHO;
  // const int UMX  = ProbClosures::UMX;
  // const int UMY  = ProbClosures::UMY;
  // const int UMZ  = ProbClosures::UMZ;
  // const int UET  = ProbClosures::UET;
  // const int face = (idir+1)*sgn;


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

  Real gradrho= sqrt(drhox*drhox+drhoy*drhoy);        

  if (nt_level > 0)
  {
    //tag cells based on density gradient
    switch (level)
    {
      case 0:
        tagfab(i,j,k) = (gradrho > 0.3);        
        break;
      case 1:
        tagfab(i,j,k) = (gradrho > 0.2);        
        break;
      default:
        tagfab(i,j,k) = (gradrho > 0.1);        
        break;
    }
  }
}
////////////////////////////////////////////////////////////////////////////////




///////////////////////////////////////////////////////////////////////////////

} // namespace PROB
#endif
