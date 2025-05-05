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

using namespace amrex;
using namespace universal_constants;

namespace PROB {

// constants  
static constexpr Real Mach     = 6.0;              // bulk Mach number
static constexpr Real Mw       = 28.96e-3;         // Molecular weight
static constexpr Real gam      = 1.4;              // Adiabatic coefficient
static constexpr Real Rgas     = gas_constant/Mw;  // gas constant
static constexpr Real Reynolds = 10000;            // bulk Reynolds number
static constexpr Real Pr       = 0.7;              // Prandtl
static constexpr Real Cv       = Rgas/(gam - 1.0);
static constexpr Real Cp       = gam*Cv;
static constexpr Real viscos   = 1.0/Reynolds;    // constant viscosity
static constexpr Real lambda   = viscos*Cp/Pr;    // constant lambda
  

//////////////////////////// Physical modelling ////////////////////////////////
struct ProbParm
{ 
  // freee-stream conditions  
  Real p_oo    = 5000.0; //[Pa] free-stream pressure (at 10 km altitude)  
  Real T_oo    = 223.0;  //[K]  free-stream temperature 
  Real rho_oo  = p_oo*Mw/(gas_constant*T_oo);
  Real c_oo    = sqrt(gam*p_oo/rho_oo);
  Real u_oo    = c_oo*Mach;  
  Real eint_oo = p_oo/ (gam - Real(1.0));
  Real kin_oo  = 0.5*rho_oo*u_oo*u_oo;
 
  // right state
  Real p_r     = 5000.0; //[Pa] free-stream pressure (at 10 km altitude)  
  Real T_r     = 223.0;  //[K]  free-stream temperature 
  Real rho_r   = p_r*Mw/(gas_constant*T_r);
  Real u_r     = 0.0;
  Real eint_r  = p_r/ (gam - Real(1.0));
 
  // centre of sphere (should be same with STL file)
  Real x0 = 1.0; Real y0 = 2.0; Real z0 = 2.0;
  // initial shock position  
  Real xshock = 0.5;  

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
  static constexpr Real C2skew=0.5,C4skew=0.0016;   // Skew symmetric default
};


// CLOSURES
typedef closures_dt<indicies_t, transport_const_t<methodparm_t>,
                    calorifically_perfect_gas_t<indicies_t>> ProbClosures;

// NUMERICAL SCHEME + EQNS TO SOLVE   (Euler/NS/Source)                 

//typedef rhs_dt<riemann_t<false, ProbClosures>, no_diffusive_t, no_source_t > ProbRHS;
typedef rhs_dt<skew_t<skewparm_t,ProbClosures>, no_diffusive_t, no_source_t > ProbRHS;
// typedef rhs_dt<skew_t<skewparm_t,ProbClosures>, viscous_t<methodparm_t, ProbClosures>, no_source_t > ProbRHS;


// IBM templates
using d_image = std::ratio<5, 5>;
typedef eib_t<1,1,d_image,ProbClosures> ProbIB;

void inline inputs() {
  
  amrex::Print() << " ****** Starting ... ******* " <<  std::endl;
  amrex::Print() << " Supersonic Flow over Sphere (IBM) " <<  std::endl;

}

//////////////////////////// Initial conditions ////////////////////////////////
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void prob_initdata (int i, int j, int k, amrex::Array4<amrex::Real> const& state,
      amrex::GeometryData const& geomdata, ProbClosures const& cls, ProbParm const& pparm) {
  
  const Real* prob_lo = geomdata.ProbLo();
  const Real* prob_hi = geomdata.ProbHi();
  const Real* dx      = geomdata.CellSize();

  Real x = prob_lo[0] + (i+0.5_rt)*dx[0];
  // Real y = prob_lo[1] + (j+0.5_rt)*dx[1];
  // Real z = prob_lo[2] + (k+0.5_rt)*dx[2];
  // local vars
  Real rhot,eint,u[3]={0.0};
  
  // initial state
  if (x < pparm.xshock) { // left of shock
    rhot =  pparm.rho_oo;
    u[0] =  pparm.u_oo;
    eint =  pparm.eint_oo;
  }
  else {                  // after shock
    rhot =  pparm.rho_r;
    u[0] =  pparm.u_r; 
    eint =  pparm.eint_r;
  }

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
  const Real x = (i+0.5_rt)*dx[0];
  const Real y = (j+0.5_rt)*dx[1];
  const Real z = (k+0.5_rt)*dx[2];
  Real xrel[3];
  // coordinate relative to object
  xrel[0]= x-pparm.x0; xrel[1]= y-pparm.y0; xrel[2]= z-pparm.z0;
  Real radius =xrel[0]*xrel[0] + xrel[1]*xrel[1] + xrel[2]*xrel[2];
  const Real Rmax = 0.5_rt*0.5_rt;const Real Rmin = 0.15_rt*0.15_rt;
  // initialize thresholds at all levels
  Real rhofluc_threshold[6] = {0.3_rt,0.6_rt,0.9_rt,1000_rt,1000_rt,1000_rt};

  // refinement first step
  if ( nt==0) {
    
    if (level==0 ) {    
      tagfab(i,j,k) = (radius < Rmax ) && (radius > Rmin);
    }

  }
  else
  {
    int URHO = ProbClosures::URHO; 
    // refine close to grads of density 
    Real drhox = std::abs(sdatafab(i+1,j,k,URHO) - sdatafab(i-1,j,k,URHO));
    Real drhoy = std::abs(sdatafab(i,j+1,k,URHO) - sdatafab(i,j-1,k,URHO));
    Real rhop  = sdatafab(i,j,k,URHO);
    Real rhofluc = std::sqrt(drhox*drhox + drhoy*drhoy)/rhop ;

    tagfab(i,j,k) = (rhofluc > rhofluc_threshold[level]);

    // always refine close to body (at all levels)
    if (ibfab(i,j,k,1)) {
      for (int ii = -1; ii <= 1; ii++) {
        for (int jj = -1; jj <= 1; jj++) {
          for (int kk = -1; kk <= 1; kk++) {
            tagfab(i+ii,j+jj,k+kk) = true;
          }
        }
      }
    }
  }

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

  switch(face)
  {
    case  -1:  // EAST      
      // inflow
      s_ext[URHO] = pparm.rho_oo;
      s_ext[UMX]  = pparm.rho_oo * pparm.u_oo;
      s_ext[UMY]  = 0.0;
      s_ext[UMZ]  = 0.0;
      s_ext[UET]  = pparm.eint_oo + pparm.kin_oo;      
      break;
    case   1:  //WEST
      break;      
    default:
      break;
  }


}

}
#endif
