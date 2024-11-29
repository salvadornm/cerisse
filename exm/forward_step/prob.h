#ifndef CNS_PROB_H_
#define CNS_PROB_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_ParmParse.H>
#include <Closures.h>
#include <RHS.h>
#include <ebm.h>
#include <walltypes.h>


// 2D Shock Reflection

using namespace amrex;

namespace PROB {

static constexpr Real Mw   = 28.96e-3;  // Molecular weight
static constexpr Real gam  = 1.4;       // Adiabatic coefficient
static constexpr Real Mach = 3.0;       // Mach Number



// problem parameters 
struct ProbParm {
  Real gam = 1.4;
  Real p_l    = 153546.24;    
  Real T_l    = 169.167;
  Real c_l    = sqrt(gam*gas_constant*T_l/Mw); 
  Real u_l    =  c_l*Mach;
  Real v_l    =  0.0;
  Real rho_l  =  p_l*Mw/(gas_constant*T_l); 
  Real eint_l =  p_l/ (gam - Real(1.0));

  Real p_r    = 153546.24;    
  Real T_r    = 153.0;
  Real c_r    = sqrt(gam*gas_constant*T_r/Mw); 
  Real u_r    = 0.0;
  Real v_r    = 0.0;
  Real rho_r  = p_r*Mw/(gas_constant*T_r); ; 
  Real eint_r = p_r/(gam - Real(1.0));;  
};

// numerical method parameters
struct methodparm_t {

  public:

  static constexpr bool dissipation = true;         // no dissipation
  static constexpr int  order = 2;                  // order numerical scheme   
  static constexpr Real C2skew=0.5,C4skew=0.0016;   // Skew symmetric default

};

inline Vector<std::string> cons_vars_names={"Xmom","Ymom","Zmom","Energy","Density"};
inline Vector<int> cons_vars_type={1,2,3,0,0};

typedef closures_dt<indicies_t, visc_suth_t, cond_suth_t,
                    calorifically_perfect_gas_t<indicies_t>> ProbClosures;

//typedef rhs_dt<rusanov_t<ProbClosures>, no_diffusive_t, no_source_t > ProbRHS;
//typedef rhs_dt<riemann_t<false, ProbClosures>, no_diffusive_t, no_source_t > ProbRHS;
typedef rhs_dt<skew_t<methodparm_t, ProbClosures>, no_diffusive_t, no_source_t > ProbRHS;

// define type of wall and EBM class

typedef adiabatic_wall_t<ProbClosures> TypeWall;

typedef ebm_t<TypeWall,ProbClosures> ProbEB;



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
  //Real y = prob_lo[1] + (j + Real(0.5)) * dx[1];
  
  // local vars
  Real rhot,eint,u[2];

  // final state
  if (x < 0.3) {
    rhot =  prob_parm.rho_l;
    u[0] =  prob_parm.u_l; u[1] =  prob_parm.v_l;
    eint =  prob_parm.eint_l;
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
      // post-shock conditions
      s_ext[URHO] = prob_parm.rho_l;
      s_ext[UMX]  = prob_parm.rho_l * prob_parm.u_l;
      s_ext[UMY]  = 0.0;
      s_ext[UMZ]  = 0.0;    
      s_ext[UET]  = prob_parm.eint_l + 0.5 *prob_parm.rho_l* prob_parm.u_l * prob_parm.u_l;
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

  Real rhot = sdatafab(i,j,k,ProbClosures::URHO);

  Real dengrad_threshold = 0.5;
  Real drhox = Math::abs(sdatafab(i+1,j,k,ProbClosures::URHO) -
   sdatafab(i,j,k,ProbClosures::URHO))/rhot;

  Real drhoy = Math::abs(sdatafab(i,j+1,k,ProbClosures::URHO) -
   sdatafab(i,j-1,k,ProbClosures::URHO))/rhot;

  Real gradrho= sqrt(drhox*drhox+drhoy*drhoy);        

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

     //tagfab(i,j,k) = (x < 1.0);

 // }
}
////////////////////////////////////////////////////////////////////////////////




///////////////////////////////////////////////////////////////////////////////

} // namespace PROB
#endif
