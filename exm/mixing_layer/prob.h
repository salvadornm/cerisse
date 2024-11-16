#ifndef CNS_PROB_H_
#define CNS_PROB_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_ParmParse.H>

#include <PelePhysics.H>
#include <ReactorBase.H>

#include <Closures.h>
#include <RHS.h>

#include <Constants.h>
#include <NumParam.h>

// Mixing Layer

using namespace amrex;


namespace PROB {


static constexpr Real Reynolds = 50.0;  // Reynolds number


// problem parameters 
struct ProbParm {

 // mid plane
 Real ymix = 0.0;

 // bottom
  Real T1 = 545;
  Real u1 = 66910;
  Real rho1;
 
  // top
  Real T2 = 1475;
  Real u2 = 115160;
  Real rho2;
 
  Real p = 942322.5;
  Real theta_w = 0.0198; // vorticity thickness
  
  // bottom
  GpuArray<Real, NUM_SPECIES> Y1 = {
      0., 0.01277243, 0.,         0., 0., 0.10136214, 0.,
      0., 0.,         0.88586543, 0., 0., 0.};  // mass fractions [-] (molefrac
                                                // 2:1:7)
  // top
   GpuArray<Real, NUM_SPECIES> Y2 = {0.,         0.01277243, 0., 0., 0.,
                                     0.10136214, 0.,         0., 0., 0.88586543,
                                     0.,         0.,         0.};

};

// numerical method parameters
struct methodparm_t {

  public:

  static constexpr int  order = 2;                  // order numerical scheme viscous
  
};

//inline Vector<std::string> cons_vars_names={"Xmom","Ymom","Zmom","Energy","Density"};
//inline Vector<int> cons_vars_type={1,2,3,0,0};


using ProbClosures = closures_dt<indicies_t, visc_suth_t, cond_suth_t,
                                 multispecies_perfect_gas_t<indicies_t>>;
using ProbRHS =
    rhs_dt<weno_t<ReconScheme::Teno5, ProbClosures>, no_diffusive_t, reactor_t<ProbClosures>>;


void inline inputs() {
  ParmParse pp;

  pp.add("cns.order_rk", 3);   // -2, 1, 2 or 3"
  pp.add("cns.stages_rk", 3);  // 1, 2 or 3

  amrex::Print() << " ****** Starting ... *******" <<  std::endl;
  amrex::Print() << " Mixing Layer  (Oct 2024)" <<  std::endl;

}

// initial condition
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
prob_initdata(int i, int j, int k, Array4<Real> const &state,
              GeometryData const &geomdata, ProbClosures const &cls,
              ProbParm const &prob_parm) {
  const Real *prob_lo = geomdata.ProbLo();
  const Real *prob_hi = geomdata.ProbHi();
  const Real *dx = geomdata.CellSize();

  //Real x = prob_lo[0] + (i + Real(0.5)) * dx[0];
  Real y = prob_lo[1] + (j + Real(0.5)) * dx[1];

  
  Real smooth= 0.5  + 0.5 * tanh(2 * y / prob_parm.theta_w);
  
  // local vars  
  Real Pt, rhot, uxt, et;
  Real Yt[NUM_SPECIES];
  
  // assign rho,P,Y based on y
  Pt = prob_parm.p;
  uxt  = prob_parm.u1*smooth   + (1.0 - smooth)*prob_parm.u2;
  rhot = prob_parm.rho1*smooth + (1.0 - smooth)*prob_parm.rho2;
  for (int n = 0; n < NUM_SPECIES; n++){
    Yt[n]   = prob_parm.Y1[n]*smooth + (1.0 - smooth)*prob_parm.Y2[n];
  }
  
  cls.RYP2E(rhot, Yt, Pt, et);

  state(i, j, k, cls.UMX) = rhot * uxt;
  state(i, j, k, cls.UMY) = Real(0.0);
  state(i, j, k, cls.UMZ) = Real(0.0);
  state(i, j, k, cls.UET) = rhot * et + Real(0.5) * rhot * uxt * uxt;
  for (int n = 0; n < NUM_SPECIES; n++) {
    state(i, j, k, cls.UFS + n) = rhot * Yt[n];
  }


  
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
  const int UFS  = ProbClosures::UFS;  
  const int face = (idir+1)*sgn;
   
  Real Pt, rhot, uxt, et;
  Real Yt[NUM_SPECIES];
  Real smooth= 0.5  + 0.5 * tanh(2 * x[1] / prob_parm.theta_w);

  switch(face)
  {
    case 1: // WEST INFLOW

      Pt = prob_parm.p;
      uxt  = prob_parm.u1*smooth   + (1.0 - smooth)*prob_parm.u2;
      rhot = prob_parm.rho1*smooth + (1.0 - smooth)*prob_parm.rho2;
      for (int n = 0; n < NUM_SPECIES; n++){
        Yt[n]   = prob_parm.Y1[n]*smooth + (1.0 - smooth)*prob_parm.Y2[n];
      }
     closures.RYP2E(rhot, Yt, Pt, et);
        
      s_ext[UMX]  = rhot*uxt;
      s_ext[UMY]  = 0.0;
      s_ext[UMZ]  = 0.0;
      s_ext[UET]  = rhot * et + Real(0.5) * rhot * uxt * uxt;
      for (int n = 0; n < NUM_SPECIES; n++) {
        s_ext[UFS + n] = rhot * Yt[n];        
      }
              
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


  Real Threshold[3] = {100.0,100.0,100.0};
  
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
