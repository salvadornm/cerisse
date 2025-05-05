#ifndef CNS_PROB_H_
#define CNS_PROB_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_ParmParse.H>
#include <Closures.h>
#include <RHS.h>
#include <ebm.h>
#include <walltypes.h>

// Hypersonic flow over a 2D Cylinder
// created by S Navarro-Martinez (2025)

using namespace amrex;

namespace PROB {

static constexpr Real Mw   = 28.96e-3;  // Molecular weight
static constexpr Real gam  = 1.4;       // Adiabatic coefficient
static constexpr Real Mach = 6.0;       // Mach Number

// problem parameters 

struct ProbParm {  
  // freestream state
  Real p_oo    = 5000.0; //[Pa] free-stream pressure (at 10 km altitude)  
  Real T_oo    = 223.0;  //[K]  free-stream temperature 
  Real rho_oo  = p_oo*Mw/(gas_constant*T_oo);
  Real c_oo    = sqrt(gam*p_oo/rho_oo);  
  Real u_oo    = c_oo*Mach;   
  Real v_oo    = 0.0;
  Real eint_oo = p_oo/ (gam - Real(1.0));  

  // right state
  Real p_r     = 5000.0; //[Pa] free-stream pressure (at 10 km altitude)  
  Real T_r     = 223.0;  //[K]  free-stream temperature 
  Real rho_r   = p_r*Mw/(gas_constant*T_r);
  Real u_r     = 0.0;   
  Real v_r     = 0.0;
  Real eint_r  = p_r/ (gam - Real(1.0));  

  // centre of cylinder (from input file)
  Real x0 = 1; Real y0 = 2;
  Real xshock = 0.75;  //[m] initial shock position
};


// numerical method parameters
struct methodparm_t {

  public:

  static constexpr bool dissipation = true;         // no dissipation
  static constexpr int  order = 4;                  // order numerical scheme   
  static constexpr Real C2skew=1.5,C4skew=0.0016;   // Skew symmetric default
  static constexpr bool solve_diffwall = false;     // solve viscous effects at walls						    
};

inline Vector<std::string> cons_vars_names={"Xmom","Ymom","Zmom","Energy","Density"};
inline Vector<int> cons_vars_type={1,2,3,0,0};

typedef closures_dt<indicies_t, visc_suth_t, cond_suth_t,
                    calorifically_perfect_gas_t<indicies_t>> ProbClosures;

// define nuemrical scheme comment/uncomment to set up 
//typedef rhs_dt<rusanov_t<ProbClosures>, no_diffusive_t, no_source_t > ProbRHS;
//typedef rhs_dt<riemann_t<false, ProbClosures>, no_diffusive_t, no_source_t > ProbRHS;
//typedef rhs_dt<skew_t<methodparm_t, ProbClosures>, no_diffusive_t, no_source_t > ProbRHS;
typedef rhs_dt<weno_t<ReconScheme::Teno5, ProbClosures>, no_diffusive_t, no_source_t > ProbRHS;


// define type of wall and EBM class
typedef adiabatic_wall_t<ProbClosures> TypeWall;
typedef ebm_t<TypeWall,methodparm_t,ProbClosures> ProbEB;


void inline inputs() {
  //	
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
  //Real y = prob_lo[1] + (j + Real(0.5)) * dx[1];
  
  // local vars
  Real rhot,eint,u[2];

  // final state
  if (x < prob_parm.xshock) {
    rhot =  prob_parm.rho_oo;
    u[0] =  prob_parm.u_oo; u[1] =  prob_parm.v_oo;
    eint =  prob_parm.eint_oo;
  }
  else {
    rhot =  prob_parm.rho_r;
    u[0] =  prob_parm.u_r; u[1] =  prob_parm.v_r;
    eint =  prob_parm.eint_r;
  }



  
  state(i, j, k, cls.URHO) = rhot;
  state(i, j, k, cls.UMX)  = rhot * u[0];
  state(i, j, k, cls.UMY)  = rhot * u[1];
  state(i, j, k, cls.UMZ)  = Real(0.0);  
  state(i, j, k, cls.UET)  = eint + Real(0.5) * rhot * (u[0] * u[0] + u[1] * u[1]);    
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
    case  2:  // SOUTH
      break;
    case  1:  // WEST
      // inflow
      s_ext[URHO] = prob_parm.rho_oo;
      s_ext[UMX]  = prob_parm.rho_oo * prob_parm.u_oo;
      s_ext[UMY]  = 0.0;
      s_ext[UMZ]  = 0.0;    
      s_ext[UET]  = prob_parm.eint_oo + 0.5 *prob_parm.rho_oo* prob_parm.u_oo * prob_parm.u_oo;
      break;
    case -1:  // EAST
      break;
    case -2:   // NORTH
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
  const Real *dx = geomdata.CellSize();
  Real x = prob_lo[0] + (i + Real(0.5)) * dx[0];
  Real y = prob_lo[1] + (j + Real(0.5)) * dx[1];
  // coordinate relative to object
  Real xrel[2];
  xrel[0]= x-prob_parm.x0; xrel[1]= y-prob_parm.y0;

  Real Rad2 = xrel[0]*xrel[0]+xrel[1]*xrel[1]; 


  Real rhot = sdatafab(i,j,k,ProbClosures::URHO);

//  Real dengrad_threshold = 0.5;
  Real drhox = Math::abs(sdatafab(i+1,j,k,ProbClosures::URHO) -
   sdatafab(i,j,k,ProbClosures::URHO))/rhot;

  Real drhoy = Math::abs(sdatafab(i,j+1,k,ProbClosures::URHO) -
   sdatafab(i,j-1,k,ProbClosures::URHO))/rhot;

  Real gradrho= sqrt(drhox*drhox+drhoy*drhoy);        

  Real Rmax = 0.3;

  //if (nt_level > 0)
 // {
    //tag cells based on density gradient
    switch (level)
    {
      case 0:
        tagfab(i,j,k) = (gradrho > 0.1);        
        break;
      case 1:
        tagfab(i,j,k) = (gradrho > 0.2);        
        break;
      default:
        tagfab(i,j,k) = (gradrho > 0.3);        
        break;
    }

    tagfab(i,j,k) = tagfab(i,j,k) || (Rad2 < Rmax*Rmax);

 // }
  
  // refine next to body


}
////////////////////////////////////////////////////////////////////////////////




///////////////////////////////////////////////////////////////////////////////

} // namespace PROB
#endif
