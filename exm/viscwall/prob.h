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
  Real dpdx = 5e-5;          // forcing term
};

// numerical method parameters
struct methodparm_t {

  public:

  static constexpr int  order = 2;              // order numerical scheme viscous
  static constexpr Real conductivity = 0.0262;  // conductivity (for constant value)
  static constexpr Real viscosity    = 0.02;  // viscosity    (for constant value)
  
};

inline Vector<std::string> cons_vars_names={"Xmom","Ymom","Zmom","Energy","Density"};
inline Vector<int> cons_vars_type={1,2,3,0,0};

typedef closures_dt<indicies_t, visc_const_t<methodparm_t>, cond_const_t<methodparm_t>,
                    calorifically_perfect_gas_t<indicies_t>> ProbClosures;

template <typename cls_t > class user_source_t;

typedef rhs_dt<no_euler_t, viscous_t<methodparm_t, ProbClosures>, user_source_t<ProbClosures> > ProbRHS;



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

  //Real x = prob_lo[0] + (i + Real(0.5)) * dx[0];
  Real y = prob_lo[1] + (j + Real(0.5)) * dx[1];

  Real yoh = y/prob_parm.L;
  
  // local vars
  Real rho,eint, T, P,u[2]; 


  // debugging
  Real const dt = 5e-3;

  // std::cout << " Info " << std::endl;

  // printf(" DIFF number = %f \n", dt*prob_parm.visc/(dx[0]*dx[0]) );
  // printf(" CFL  number = %f \n", dt*prob_parm.u/dx[0] );
  
  // exit(0);

  //

  
  // initial conditioms  
  P    = prob_parm.p; rho = prob_parm.rho;
  //u[0] = prob_parm.u; 
  u[1] = Real(0.0);
  u[0] = prob_parm.u*(1 - yoh)*yoh;

  // T and E from P and T
  T =   P/(cls.Rspec*rho);
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

  const int URHO = ProbClosures::URHO;
  const int UMX  = ProbClosures::UMX;
  const int UMY  = ProbClosures::UMY;
  const int UMZ  = ProbClosures::UMZ;
  const int UET  = ProbClosures::UET;
  const int face = (idir+1)*sgn;

  switch(face)
  {
    case 2: // SOUTH adiabatic wall
      s_ext[URHO] = s_int[URHO];
      s_ext[UMX]  = -s_int[UMX];
      s_ext[UMY]  = -s_int[UMY];
      s_ext[UMZ]  = -s_int[UMZ];
      s_ext[UET]  = s_int[UET];
      break;
    case -2: // NORTH adiabatic wall
      s_ext[URHO] = s_int[URHO];
      s_ext[UMX]  = -s_int[UMX];
      s_ext[UMY]  = -s_int[UMY];
      s_ext[UMZ]  = -s_int[UMZ];
      s_ext[UET]  = s_int[UET];
    default:
      break;  
  }

}
///////////////////////////////SOURCE TERM /////////////////////////////////////
template <typename cls_t>
class user_source_t {
  public:
  void inline src(const amrex::MFIter &mfi,
                  const amrex::Array4<const amrex::Real> &prims,
                  const amrex::Array4<amrex::Real> &rhs, const cls_t *cls_d,
                  amrex::Real dt){

    const Box bx = mfi.tilebox();

    ProbParm const prob_parm;


    amrex::ParallelFor(bx,
      [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
      const auto& cls = *cls_d;

      //  Pressure Source           
      //Real rho = prims(i, j, k, cls.QRHO);
      rhs(i,j,k,cls.UMX) += prob_parm.dpdx;      
      });

  };
};

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
