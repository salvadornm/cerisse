#ifndef CNS_PROB_H_
#define CNS_PROB_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_ParmParse.H>
#include <Closures.h>
#include <RHS.h>

#include <Constants.h>
#include <NumParam.h>

// 2D periodic channel laminar flow
// Pure diffusion case taken from HAMISH
// https://www.ukctrf.com/index.php/benchmarking-of-the-new-software/

using namespace amrex;

namespace PROB {

// problem parameters 
struct ProbParm {
  Real Re   = 50.0;          // Reynodls number
  Real L    =   1.0;         // channel height    
  Real p    =  pres_atm2si; 
  Real rho  = 1.0;
  Real u    = 1.0;
  Real visc = rho*u*L/Re;    // default  0.02
};

// numerical method parameters
struct methodparm_t {

  public:

  static constexpr int  order = 2;              // order numerical scheme viscous
  static constexpr Real conductivity = 0.0262;  // conductivity (for constant value)
  static constexpr Real viscosity   = 1.85e-5;  // viscosity    (for constant value)
  
};

inline Vector<std::string> cons_vars_names={"Xmom","Ymom","Zmom","Energy","Density"};
inline Vector<int> cons_vars_type={1,2,3,0,0};

typedef closures_dt<indicies_t, visc_const_t<methodparm_t>, cond_const_t<methodparm_t>,
                    calorifically_perfect_gas_t<indicies_t>> ProbClosures;

typedef rhs_dt<no_euler_t, viscous_t<methodparm_t, ProbClosures>, no_source_t > ProbRHS;



void inline inputs() {
  ParmParse pp;

  pp.add("cns.order_rk", 3);   // -2, 1, 2 or 3"
  pp.add("cns.stages_rk", 3);  // 1, 2 or 3

  amrex::Print() << " ****** Starting ... *******" <<  std::endl;
  amrex::Print() << " 2D Periodic Channel Flow  (Oct 2024)" <<  std::endl;

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
  P    = prob_parm.p; rho = prob_parm.rho;
  u[0] = prob_parm.u; u[1] = Real(0.0);

  // T and E from P and T
  T =   P/(cls.Rspec*rho);
  eint =  P / (cls.gamma - Real(1.0));

  state(i, j, k, cls.URHO) = rho;
  state(i, j, k, cls.UMX)  = rho * u[0];
  state(i, j, k, cls.UMY)  = rho * u[1];
  state(i, j, k, cls.UMZ)  = Real(0.0);  
  state(i, j, k, cls.UET)  = eint + Real(0.5) * rho * (u[0] * u[0] + u[1] * u[1]);    
}

/////////////////////////////// Source /////////////////////////////////////////////


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

  // south adiabatic wall

  // north adiabatic wall


}
///////////////////////////////AMR//////////////////////////////////////////////
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
user_tagging(int i, int j, int k, int nt_level, auto &tagfab,
             const auto &sdatafab, const auto &geomdata,
             const ProbParm &prob_parm, int level) {


  Real Threshold[3] = {0.1,0.05,0.02};
  
  Real rhop = sdatafab(i,j,k,ProbClosures::URHO);
  Real drhox = Math::abs(sdatafab(i+1,j,k,ProbClosures::URHO) -
              sdatafab(i-1,j,k,ProbClosures::URHO));
  Real drhoy = Math::abs(sdatafab(i,j+1,k,ProbClosures::URHO) -
              sdatafab(i,j-1,k,ProbClosures::URHO));
  Real gradrho= sqrt(drhox*drhox+drhoy*drhoy)/rhop;

  if (nt_level> 0)
  {
    //tag cells based on density  values       
    tagfab(i,j,k) = (gradrho > Threshold[level]);        
      
  }
}
////////////////////////////////////////////////////////////////////////////////




///////////////////////////////////////////////////////////////////////////////

} // namespace PROB
#endif
