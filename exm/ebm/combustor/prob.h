#ifndef CNS_PROB_H_
#define CNS_PROB_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_ParmParse.H>
#include <Closures.h>
#include <RHS.h>

#if CNS_USE_EB    
#include <ebm.h>
#include <walltypes.h>
#endif

// NTU Combustor-type for demonstration purposes
// created by S Dupre and S Navarro-Martinez (2025)


using namespace amrex;

namespace PROB {

static constexpr Real Mw   = 28.96e-3;  // Molecular weight
static constexpr Real gam  = 1.4;       // Adiabatic coefficient

// problem parameters 

struct ProbParm {  
  // inflow state
  Real p_inflow    = pres_atm2si; //[Pa] inflow pressure (1 atm)  
  Real T_inflow    = 298.0;  //[K]  
  Real rho_inflow  = p_inflow*Mw/(gas_constant*T_inflow);
  Real vel_in[3]= {0.0,0.0,10.0}; // array of inflow velocity [m/s]
  Real eint_inflow = p_inflow/ (gam - Real(1.0));  

  Real rhou_inflow = rho_inflow*vel_in[0];
  Real rhov_inflow = rho_inflow*vel_in[1];
  Real rhow_inflow = rho_inflow*vel_in[2];
  Real kin = Real(0.5) * rho_inflow *
      (vel_in[0]*vel_in[0]+vel_in[1]*vel_in[1]+vel_in[2]*vel_in[2]);
  Real rhoe_inflow = eint_inflow + kin;  
  
  // inside combustor state
  Real p_0     = pres_atm2si; //[Pa] inflow pressure (1 atm) 
  Real T_0     = 300.0;  //[K]  
  Real rho_0   = p_0*Mw/(gas_constant*T_0);
  Real vel_0[3]= {0.0,0.0,0.0}; // array of inside velocity [m/s]
  Real eint_0  = p_0/ (gam - Real(1.0));  

  // geometry auxiliary
  Real zin = 0.01;
};


// numerical method parameters
struct methodparm_t {

  public:

  static constexpr bool dissipation = true;         // no dissipation
  static constexpr int  order = 4;                  // order numerical scheme   
  static constexpr Real C2skew=0.5,C4skew=0.0016;   // Skew symmetric default

};

inline Vector<std::string> cons_vars_names={"Xmom","Ymom","Zmom","Energy","Density"};
inline Vector<int> cons_vars_type={1,2,3,0,0};

typedef closures_dt<indicies_t, visc_suth_t, cond_suth_t,
                    calorifically_perfect_gas_t<indicies_t>> ProbClosures;

// define nuemrical scheme comment/uncomment to set up 
//typedef rhs_dt<rusanov_t<ProbClosures>, no_diffusive_t, no_source_t > ProbRHS;
//typedef rhs_dt<riemann_t<false, ProbClosures>, no_diffusive_t, no_source_t > ProbRHS;
typedef rhs_dt<skew_t<methodparm_t, ProbClosures>, no_diffusive_t, no_source_t > ProbRHS;
//typedef rhs_dt<weno_t<ReconScheme::Teno5, ProbClosures>, no_diffusive_t, no_source_t > ProbRHS;


// define type of wall and EBM class
#if CNS_USE_EB    
typedef adiabatic_wall_t<ProbClosures> TypeWall;
typedef ebm_t<TypeWall,ProbClosures> ProbEB;
#endif

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
  Real y = prob_lo[1] + (j + Real(0.5)) * dx[1];
  Real z = prob_lo[2] + (k + Real(0.5)) * dx[2];

  Real rad = sqrt(x*x + y*y);

  // local vars
  Real rhot,eint,u[3];

  if (z < prob_parm.zin) {
    rhot =  prob_parm.rho_inflow;    
    for(int idim=0;idim < AMREX_SPACEDIM;idim++) {u[idim]=prob_parm.vel_in[idim];}
    eint =  prob_parm.eint_inflow;
  }
  else {
    rhot =  prob_parm.rho_0;    
    for(int idim=0;idim < AMREX_SPACEDIM;idim++) {u[idim]=prob_parm.vel_0[idim];}
    eint =  prob_parm.eint_0;
  }
  
  Real kin = Real(0.5) * rhot * (u[0] * u[0] + u[1] * u[1] + u[2]*u[2]);
  state(i, j, k, cls.URHO) = rhot;
  state(i, j, k, cls.UMX)  = rhot * u[0];
  state(i, j, k, cls.UMY)  = rhot * u[1];
  state(i, j, k, cls.UMZ)  = rhot * u[2];  
  state(i, j, k, cls.UET)  = eint + kin;    
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
   
  const int face = (idir+1)*sgn; // +/-1 (1D) +/- 2 (2D) +/- 3 (3D)

  switch(face)
  {
    case  3:  // LEFT
      // inflow
      s_ext[URHO] = prob_parm.rho_inflow;
      s_ext[UMX]  = prob_parm.rhou_inflow;
      s_ext[UMY]  = prob_parm.rhov_inflow;
      s_ext[UMZ]  = prob_parm.rhow_inflow;
      s_ext[UET]  = prob_parm.rhoe_inflow;
    case  2:  // SOUTH
      break;
    case  1:  // WEST
      break;
    case -1:  // EAST
      break;
    case -2:  // NORTH
      break;
    case -3:   //RIGHT 
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
  Real z = prob_lo[2] + (k + Real(0.5)) * dx[2];

  const int URHO= ProbClosures::URHO;

  // compoute | d rho | normalised with rho
  Real o_over_rhot = Real(1.0)/sdatafab(i,j,k,URHO);

  Real drhox = Math::abs(sdatafab(i+1,j,k,URHO) - sdatafab(i-1,j,k,URHO));
  Real drhoy = Math::abs(sdatafab(i,j+1,k,URHO) - sdatafab(i,j-1,k,URHO));
  Real drhoz = Math::abs(sdatafab(i,j,k+1,URHO) - sdatafab(i,j,k-1,URHO));

  Real gradrho= Real(0.5)*sqrt(drhox*drhox+drhoy*drhoy)*o_over_rhot;        


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
    
  // refine next to body


}
///////////////////////////////////////////////////////////////////////////////

} // namespace PROB

#endif
