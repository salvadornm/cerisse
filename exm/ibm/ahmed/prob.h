#ifndef CNS_PROB_H_
#define CNS_PROB_H_

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_ParmParse.H>
#include <AMReX_AmrLevel.H>
#include <Closures.h>
#include <RHS.h>
#include <eib.h>
#include <ratio>
#include <Constants.h>
#include <ib_walltypes.h>


using namespace amrex;
using namespace universal_constants;

namespace PROB {

// constants  
static constexpr Real Mw       = 28.96e-3;         // Molecular weight
static constexpr Real gam      = 1.4;              // Adiabatic coefficient
static constexpr Real Rgas     = gas_constant/Mw;  // gas constant
static constexpr Real Reynolds = 10000;            // bulk Reynolds number (approx)
static constexpr Real Pr       = 0.7;              // Prandtl
static constexpr Real Cv       = Rgas/(gam - 1.0);
static constexpr Real Cp       = gam*Cv;
static constexpr Real viscos   = 1.0/Reynolds;    // constant viscosity
// static constexpr Real viscos   = 1.85e-5;    // constant viscosity
static constexpr Real lambda   = viscos*Cp/Pr;    // constant lambda
  

//////////////////////////// Physical modelling ////////////////////////////////
struct ProbParm
{ 
  // freee-stream conditions  
  Real p_oo    = 100000.0; //[Pa] free-stream pressure   
  Real T_oo    = 300.0;    //[K]  free-stream temperature 
  Real u_oo    = 40.0;     // m/s
  Real rho_oo  = p_oo/(Rgas*T_oo);  
  Real eint_oo = rho_oo*Cv*T_oo;
  Real kin_oo  = 0.5*rho_oo*u_oo*u_oo;
 
  // box that contains body
  Real xbox_min = -0.2, xbox_max = 1.2;
  Real ybox_min = -0.4, ybox_max = 0.4;
  Real zbox_min = 0,    zbox_max = 0.4;

  Real x0=0,y0=0,z0=0; // vehicle position 
  Real bl= 0.06; // 60 mm BL

};

//  parameters for viscous solver and conductivity/viscosity
struct methodparm_t {

  public:

  static constexpr int  order = 2;                  // order numerical scheme   
  static constexpr Real conductivity = lambda;       // conductivity (for constant value)
  static constexpr Real viscosity    = viscos;       // viscosity    (for constant value)
};


// parameters for skew-symmetric method
struct skewparm_t {

  public:

  static constexpr bool dissipation = true;         // no dissipation
  static constexpr int  order = 4;                  // order numerical scheme   
  static constexpr Real C2skew=0.5,C4skew=0.016;   // Skew symmetric values 
};

struct ibmparm_t {

  public:

  static constexpr int  interp_order = 1;
  static constexpr int  extrap_order = 1;
  static constexpr Real alpha= 0.6;      
};


// CLOSURES
typedef closures_dt<indicies_t, transport_const_t<methodparm_t>,
                    calorifically_perfect_gas_t<indicies_t>> ProbClosures;

// NUMERICAL SCHEME + EQNS TO SOLVE   (Euler/NS/Source)                 

//typedef rhs_dt<riemann_t<false, ProbClosures>, no_diffusive_t, no_source_t > ProbRHS;
//typedef rhs_dt<skew_t<skewparm_t,ProbClosures>, no_diffusive_t, no_source_t > ProbRHS;
//typedef rhs_dt<skew_t<skewparm_t,ProbClosures>, viscous_t<methodparm_t, ProbClosures>, no_source_t > ProbRHS;
typedef rhs_dt<riemann_t<false, ProbClosures>, viscous_t<methodparm_t, ProbClosures>, no_source_t > ProbRHS;


// IBM templates

typedef ibm_adiabatic_noslip_wall_t<ibmparm_t,ProbClosures> TypeWall;
typedef eib_t<TypeWall,ibmparm_t,ProbClosures> ProbIB;


void inline inputs() {
  
  amrex::Print() << " ****** Starting ... ******* " <<  std::endl;
  amrex::Print() << "  Flow over Ahmed vehicle (IBM) " <<  std::endl;

}

//////////////////////////// Initial conditions ////////////////////////////////
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void prob_initdata (int i, int j, int k, amrex::Array4<amrex::Real> const& state,
      amrex::GeometryData const& geomdata, ProbClosures const& cls, ProbParm const& pparm) {
  
  const Real* prob_lo = geomdata.ProbLo();
  // const Real* prob_hi = geomdata.ProbHi(); //
  const Real* dx      = geomdata.CellSize();

  Real x = prob_lo[0] + (i+0.5_rt)*dx[0];
  // Real y = prob_lo[1] + (j+0.5_rt)*dx[1];
  // Real z = prob_lo[2] + (k+0.5_rt)*dx[2];
  // local vars
  Real rhot,eint,u[3]={0.0};
  
  // initial state
  rhot =  pparm.rho_oo;
  u[0] =  pparm.u_oo;
  eint =  pparm.eint_oo;

 // u[0]=0.0;
 
  state(i, j, k, cls.URHO) = rhot;
  state(i, j, k, cls.UMX)  = rhot * u[0];
  state(i, j, k, cls.UMY)  = Real(0.0);
  state(i, j, k, cls.UMZ)  = Real(0.0);
  state(i, j, k, cls.UET)  = eint + Real(0.5) * rhot * u[0] * u[0] ; 

}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE 
void user_tagging(int i, int j, int k, int nt, auto& tagfab, const auto &sdatafab, 
                  const Array4<bool>&ibfab, const auto& geomdata, 
                  const ProbParm& pparm , int level) {

  const Real* dx  = geomdata.CellSize();
  const Real* prob_lo = geomdata.ProbLo();
  const Real x = prob_lo[0] + (i+0.5_rt)*dx[0];
  const Real y = prob_lo[1] + (j+0.5_rt)*dx[1];
  const Real z = prob_lo[2] + (k+0.5_rt)*dx[2];
  Real xrel[3] = {x-pparm.x0,y-pparm.y0,z-pparm.z0};
  // coordinate relative to object
  //xrel[0]= x-pparm.x0; xrel[1]= y-pparm.y0; xrel[2]= z-pparm.z0;

  // manual refinement surroudings of car 
    
   // if (level==0 ) {    
      tagfab(i,j,k) =   (xrel[0] < pparm.xbox_max) && (xrel[0] > pparm.xbox_min) &&
                        (xrel[1] < pparm.ybox_max) && (xrel[1] > pparm.ybox_min) &&
                        (xrel[2] < pparm.zbox_max) ;
    //}  

}
//////////////////////////// Boundary conditions ///////////////////////////////
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
bcnormal(const amrex::Real x[AMREX_SPACEDIM], amrex::Real dratio, const amrex::Real s_int[ProbClosures::NCONS],
         const amrex::Real s_refl[ProbClosures::NCONS], amrex::Real s_ext[ProbClosures::NCONS],
         const int idir, const int sgn, const amrex::Real time,
         amrex::GeometryData const& /*geomdata*/,  ProbClosures const& closures, ProbParm const& pparm)  
{
  const int URHO = ProbClosures::URHO;
  const int UMX  = ProbClosures::UMX;
  const int UMY  = ProbClosures::UMY;
  const int UMZ  = ProbClosures::UMZ;
  const int UET  = ProbClosures::UET;
  const int face = (idir+1)*sgn;
  

  const int z = x[3]; //vertical coordinate 


  Real u_inflow = pparm.u_oo;

  // Real zz = z/pparm.bl;

  // if (zz < 1)
  // {
  //   u_inflow = pparm.u_oo*zz*zz;
  // } 

  Real kin =  0.5*pparm.rho_oo * u_inflow * u_inflow;

  switch(face)
  {
    case  -1:  // EAST      
      break;
    case  1:  // WEST
      // inflow
      s_ext[URHO] = pparm.rho_oo;
      s_ext[UMX]  = pparm.rho_oo * u_inflow;
      s_ext[UMY]  = 0.0;
      s_ext[UMZ]  = 0.0;
      s_ext[UET]  = pparm.eint_oo + kin;      
      break;  
    default:
      break;
  }


}

}
#endif
