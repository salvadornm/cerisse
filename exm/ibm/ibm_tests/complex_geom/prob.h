#ifndef CNS_PROB_H_
#define CNS_PROB_H_

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_ParmParse.H>
#include <AMReX_AmrLevel.H>
#include <Closures.h>
#include <RHS.h>
#include <ibm_solver.h>
#include <ratio>
#include <Constants.h>
#include <ibm_walltypes.h>


using namespace amrex;
using namespace universal_constants;

namespace PROB {

// constants  
static constexpr Real Mach     = 4.0;              // bulk Mach number
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
  Real p_oo    = 100000.0; //[Pa] free-stream pressure   
  Real T_oo    = 300.0;    //[K]  free-stream temperature 
  Real rho_oo  = p_oo/(Rgas*T_oo);
  Real c_oo    = sqrt(gam*Rgas*T_oo);
  Real u_oo    = c_oo*Mach;  
  Real eint_oo = rho_oo*Cv*T_oo;
  Real kin_oo  = 0.5*rho_oo*u_oo*u_oo;
 
  // right state
  Real p_r     = p_oo; //[Pa] free-stream pressure (at 10 km altitude)  
  Real T_r     = T_oo;  //[K]  free-stream temperature 
  Real rho_r   = p_r/(Rgas*T_r);
  Real u_r     = 0.0;
  Real eint_r  = rho_r*Cv*T_r;
 
  // centre of sphere (should be same with STL file)
  Real x0 = 0.0; Real y0 = 0.0; 
#if (AMREX_SPACEDIM == 3)
  Real z0 = 0.0;
#endif
  // initial shock position  
  Real xshock = -0.9;  

};

//  parameters for viscous solver and conductivity/viscosity
struct methodparm_t {

  public:

  static constexpr int  order = 2;                   // order numerical scheme   
  static constexpr Real conductivity = lambda;       // conductivity (for constant value)
  static constexpr Real viscosity    = viscos;       // viscosity    (for constant value)
  static constexpr bool use_LES = false;             // LES model switch
};


// parameters for skew-symmetric method
struct skewparm_t {

  public:

  static constexpr bool dissipation = true;         // no dissipation
  static constexpr int  order = 4;                  // order numerical scheme   
  static constexpr Real C2skew=1.5,C4skew=0.016;   // Skew symmetric values 
};

struct ibmparm_t {

  public:

  static constexpr int  interp_order = 1;
  static constexpr int  extrap_order = 1;
  static constexpr Real alpha = 0.6;   
  
  static constexpr int  interp_order_surf = 1;
  static constexpr int  extrap_order_surf = 1;
  static constexpr Real alpha_surf = 0.6; 

  static constexpr int  ghost_layers = 1;
  static constexpr bool interior_is_solid = true; // true: interior of geometry is solid; false: interior is fluid
};



// CLOSURES
typedef closures_dt<indicies_t, transport_const_t<methodparm_t>,
                    calorifically_perfect_gas_t<indicies_t>> ProbClosures;

// NUMERICAL SCHEME + EQNS TO SOLVE   (Euler/NS/Source)                 

//typedef rhs_dt<rusanov_t<ProbClosures>, no_diffusive_t, no_source_t > ProbRHS;
typedef rhs_dt<riemann_t<false, ProbClosures>, no_diffusive_t, no_source_t > ProbRHS;
//typedef rhs_dt<skew_t<skewparm_t,ProbClosures>, no_diffusive_t, no_source_t > ProbRHS;
//typedef rhs_dt<skew_t<skewparm_t,ProbClosures>, viscous_t<methodparm_t, ProbClosures>, no_source_t > ProbRHS;


// IBM templates
//using d_image = std::ratio<5, 5>;
//typedef ibm_solver_t<1,1,d_image,ProbClosures> ProbIB;


typedef ibm_adiabatic_noslip_wall_t<ibmparm_t,ProbClosures> TypeWall;
typedef ibm_solver_t<TypeWall,ibmparm_t,ProbClosures> ProbIB;

// Static geometry: no update needed
inline void update_geometry(Real /*time*/,
                            Vector<GeomType>& /*geom_a*/,
                            int /*ngeom*/) {}

void inline inputs() {
  
  amrex::Print() << " ****** Starting ... ******* " <<  std::endl;
  amrex::Print() << " Supersonic Flow over Sphere (IBM) " <<  std::endl;

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
#if (AMREX_SPACEDIM == 3)
  state(i, j, k, cls.UMZ)  = Real(0.0);
#endif
  state(i, j, k, cls.UET)  = eint + Real(0.5) * rhot * u[0] * u[0] ; 

}



AMREX_GPU_DEVICE AMREX_FORCE_INLINE 
void user_tagging(int i, int j, int k, int nt, auto& tagfab, const auto &sdatafab, 
                  const Array4<uint8_t>&ibfab, const auto& geomdata, 
                  const ProbParm& pparm , int level) {
  
  const Real* prob_lo = geomdata.ProbLo();                  
  const Real* dx  = geomdata.CellSize();
  const Real x =  prob_lo[0] + (i+0.5_rt)*dx[0];
  const Real y =  prob_lo[1] + (j+0.5_rt)*dx[1];
#if (AMREX_SPACEDIM == 3)
  const Real z =  prob_lo[2] + (k+0.5_rt)*dx[2];
  // coordinate relative to object
  Real xrel[3];  
  xrel[0]= x-pparm.x0; xrel[1]= y-pparm.y0; xrel[2]= z-pparm.z0;
  Real radius = sqrt(xrel[0]*xrel[0] + xrel[1]*xrel[1] + xrel[2]*xrel[2]);
#else
  Real xrel[2];
  xrel[0]= x-pparm.x0; xrel[1]= y-pparm.y0;
  Real radius = sqrt(xrel[0]*xrel[0] + xrel[1]*xrel[1]);
#endif
  const Real Rmax = 0.5_rt;const Real Rmin = 0.2_rt;
  // initialize thresholds at all levels
  Real rhofluc_threshold[6] = {0.3_rt,0.6_rt,0.9_rt,1000_rt,1000_rt,1000_rt};


  //tagfab(i,j,k) = (radius < Rmax ) && (radius > Rmin);
  
  // refine close to grads of density 
  int URHO = ProbClosures::URHO; 
  Real rhop  = sdatafab(i,j,k,URHO);
  Real drhox = std::abs(sdatafab(i+1,j,k,URHO) - sdatafab(i-1,j,k,URHO));
  Real drhoy = std::abs(sdatafab(i,j+1,k,URHO) - sdatafab(i,j-1,k,URHO));
#if (AMREX_SPACEDIM == 3)
  Real drhoz = std::abs(sdatafab(i,j,k+1,URHO) - sdatafab(i,j,k-1,URHO));  
  Real rhofluc = std::sqrt(drhox*drhox + drhoy*drhoy  + drhoz*drhoz)/rhop ;
#else
  Real rhofluc = std::sqrt(drhox*drhox + drhoy*drhoy)/rhop ;
#endif

  tagfab(i,j,k) = (rhofluc > rhofluc_threshold[level]);

   // refine close to body (at all levels)
    if (ibfab(i,j,k,1)) {
      for (int ii = -1; ii <= 1; ii++) {
        for (int jj = -1; jj <= 1; jj++) {
#if (AMREX_SPACEDIM == 3)
          for (int kk = -1; kk <= 1; kk++) {
            tagfab(i+ii,j+jj,k+kk) = true;
          }
#else
          tagfab(i+ii,j+jj,k) = true;
#endif
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
    case  -1:  // x-hi (right) — outflow
      break;
    case   1:  // x-lo (left) — inflow
      s_ext[URHO] = pparm.rho_oo;
      s_ext[UMX]  = pparm.rho_oo * pparm.u_oo;
      s_ext[UMY]  = 0.0;
#if (AMREX_SPACEDIM == 3)
      s_ext[UMZ]  = 0.0;
#endif
      s_ext[UET]  = pparm.eint_oo + pparm.kin_oo;      
      break;  
    default:
      break;
  }


}

}
#endif
