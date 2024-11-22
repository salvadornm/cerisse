#ifndef CNS_PROB_H_
#define CNS_PROB_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_ParmParse.H>

#include <Closures.h>
#include <RHS.h>

#include <Constants.h>
#include <NumParam.h>

// Mixing Layer

using namespace amrex;


namespace PROB {

static constexpr Real Mw  = 28.96e-3;  // Molecular weight
static constexpr Real gam = 1.4;       // Adiabatic coefficient


// problem parameters 
struct ProbParm {

  // mid plane
  Real ymix = 0.0;
  
  // pressure
  Real p = 94232.25;

  // top
  Real T1 = 545;
  Real u1 = 669.10;
  Real rho1 = p*Mw/(gas_constant*T1);
  Real c1 = sqrt(gam*p/rho1);
 
  // bottom
  Real T2 = 1475;
  Real u2 = 1151.6;
  Real rho2 = p*Mw/(gas_constant*T2);
  Real c2 = sqrt(gam*p/rho2);

  Real uc = (u2*c2+u1*c1)/(c1+c2);
 
 
  Real theta_w = 1.98e-4; // vorticity thickness  ~ 0.2 mm
  
  //Real theta_w = 0.2;
  
};

inline Vector<std::string> cons_vars_names={"Xmom","Ymom","Zmom","Energy","Density"};
inline Vector<int> cons_vars_type={1,2,3,0,0};


using ProbClosures = closures_dt<indicies_stat_t, visc_suth_t, cond_suth_t,
                                 calorifically_perfect_gas_t<indicies_t>>;
using ProbRHS =
    rhs_dt<weno_t<ReconScheme::Teno5, ProbClosures>, viscous_t< defaultparm_t, ProbClosures> , no_source_t >;


void inline inputs() {
  ParmParse pp;

  amrex::Print() << " ****** Starting ... *******" <<  std::endl;
  amrex::Print() << " Mixing Layer  (Oct 2024)" <<  std::endl;
  amrex::Print() << " non-reactive" <<  std::endl;
  amrex::Print() << " ***************************" <<  std::endl; 
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
  Real r = y*y/(prob_parm.theta_w*prob_parm.theta_w);

  // at top (smooth =1)  at bootomo (smooth = 2)
  Real smooth= 0.5  + 0.5 * tanh(2 * y / prob_parm.theta_w);



  
  // local vars  
  Real Pt, rhot, uxt, et;
  
  // assign rho,P,Y based on y
  Pt = prob_parm.p;
  uxt  = prob_parm.u1*smooth   + (1.0 - smooth)*prob_parm.u2;
  rhot = prob_parm.rho1*smooth + (1.0 - smooth)*prob_parm.rho2;
  
  // fluctuation in y, added close to cemtre
  Real uyt = 0.1*prob_parm.uc*sin(8*2.0*M_PI*x/prob_hi[0])*exp(-r);


  // rho e
  et =  Pt / (cls.gamma - Real(1.0));  
   
  
  Real kin = 0.5*rhot*(uxt*uxt+uyt*uyt);



  state(i, j, k, cls.URHO) = rhot;
  state(i, j, k, cls.UMX)  = rhot * uxt;
  state(i, j, k, cls.UMY)  = rhot * uyt;
  state(i, j, k, cls.UMZ)  = Real(0.0);
  state(i, j, k, cls.UET)  = et + kin;

  
}

/////////////////////////////// BC /////////////////////////////////////////////
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
bcnormal(const Real x[AMREX_SPACEDIM], Real dratio, const Real s_int[ProbClosures::NCONS],
         const Real s_refl[ProbClosures::NCONS], Real s_ext[ProbClosures::NCONS], const int idir,
         const int sgn, const Real time, GeometryData const & /*geomdata*/,
         ProbClosures const &closures, ProbParm const &prob_parm) {

  const int URHO = ProbClosures::URHO;
  const int UMX  = ProbClosures::UMX;
  const int UMY  = ProbClosures::UMY;
  const int UMZ  = ProbClosures::UMZ;
  const int UET  = ProbClosures::UET;
  const int face = (idir+1)*sgn;
   
  Real Pt, rhot, uxt, et;
  Real smooth= 0.5  + 0.5 * tanh(2 * x[1] / prob_parm.theta_w);

  switch(face)
  {
    case 1: // WEST INFLOW

      Pt = prob_parm.p;
      uxt  = prob_parm.u1*smooth   + (1.0 - smooth)*prob_parm.u2;
      rhot = prob_parm.rho1*smooth + (1.0 - smooth)*prob_parm.rho2;
      et =  Pt / (closures.gamma - Real(1.0));  

      s_ext[URHO] = rhot;  
      s_ext[UMX]  = rhot*uxt;
      s_ext[UMY]  = 0.0;
      s_ext[UMZ]  = 0.0;
      s_ext[UET]  = et + Real(0.5) * rhot * uxt * uxt;
                    
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


  const Real *prob_lo = geomdata.ProbLo();
  const Real *prob_hi = geomdata.ProbHi();
  const Real *dx = geomdata.CellSize();

  Real y = prob_lo[1] + (j + Real(0.5)) * dx[1];

  //Real Threshold[3]  = {0.05,0.1,0.1};
  Real YThreshold[3] = {20*prob_parm.theta_w,5*prob_parm.theta_w,2*prob_parm.theta_w};
  
  // Real rhop = sdatafab(i,j,k,ProbClosures::URHO);
  // Real drhox = Math::abs(sdatafab(i+1,j,k,ProbClosures::URHO) -
  //             sdatafab(i-1,j,k,ProbClosures::URHO));
  // Real drhoy = Math::abs(sdatafab(i,j+1,k,ProbClosures::URHO) -
  //             sdatafab(i,j-1,k,ProbClosures::URHO));
  // Real gradrho= sqrt(drhox*drhox+drhoy*drhoy)/rhop;

//  if (nt_level> 0)
//  {
    //tag cells based on density  values       
    // tagfab(i,j,k) = (gradrho > Threshold[level]);

    tagfab(i,j,k) = (Math::abs(y) < YThreshold[level]);            
      
 // }
}
////////////////////////////////////////////////////////////////////////////////




///////////////////////////////////////////////////////////////////////////////

} // namespace PROB
#endif
