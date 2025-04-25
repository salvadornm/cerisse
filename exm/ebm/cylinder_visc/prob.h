#ifndef CNS_PROB_H_
#define CNS_PROB_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_ParmParse.H>
#include <Closures.h>
#include <RHS.h>
#include <ebm.h>
#include <walltypes.h>

// Viscous flow over a 2D square Cylinder It will create large Von-Karman vortex
// created by S Navarro-Martinez (2025)

using namespace amrex;

namespace PROB {

static constexpr Real Mw   = 28.96e-3;    // Molecular weight
static constexpr Real gam  = 1.4;         // Adiabatic coefficient
static constexpr Real Mach = 0.5;         // Mach Number
static constexpr Real Reynolds = 100.0;   // Reynolds number


// problem parameters 

struct ProbParm {  
  // freestream state
  Real p_oo    = 153.0; //[Pa] free-stream pressure (at 10 km altitude)  
  Real T_oo    = 57.8;  //[K]  free-stream temperature 
  Real rho_oo  = p_oo*Mw/(gas_constant*T_oo);
  Real c_oo    = sqrt(gam*p_oo/rho_oo);  
  Real u_oo    = c_oo*Mach;   
  Real v_oo    = 0.0;
  Real eint_oo = p_oo/ (gam - Real(1.0));  

  // right state
  Real p_r     = p_oo; //[Pa] free-stream pressure (at 10 km altitude)  
  Real T_r     = T_oo;  //[K]  free-stream temperature 
  Real rho_r   = p_r*Mw/(gas_constant*T_r);
  Real u_r     = 0.0;   
  Real v_r     = 0.0;
  Real eint_r  = p_r/ (gam - Real(1.0));  

  // centre of object
  Real x0 = 20;
  Real y0 = 20;
};


// numerical method parameters
struct num_method_param {

  public:

  static constexpr bool dissipation = true;         // no dissipation
  static constexpr int  order = 4;                  // order numerical scheme   
  static constexpr Real C2skew=0.5,C4skew=0.016;   // Skew symmetric default

};

// viscosisty parameters
struct visc_model_param {

  public:

  static constexpr int  order = 4;                   // order numerical scheme viscous
  static constexpr Real conductivity = 0.03;         // conductivity (constant)
  static constexpr Real viscosity    = 1.0/Reynolds; // viscosity    (constant)

};

struct wall_param {

  public:
  
  static constexpr Real Twall = 300;                // wall temperature 
  static constexpr bool solve_diffwall = true;      // solve viscous effects at walls

};



inline Vector<std::string> cons_vars_names={"Xmom","Ymom","Zmom","Energy","Density"};
inline Vector<int> cons_vars_type={1,2,3,0,0};

// define thermodynamic closure
typedef closures_dt<indicies_t, transport_const_t<visc_model_param>, calorifically_perfect_gas_t<indicies_t> > ProbClosures;

// define numerical scheme  
typedef rhs_dt<skew_t<num_method_param, ProbClosures>, viscous_t<visc_model_param, ProbClosures>, no_source_t > ProbRHS;
//typedef rhs_dt<skew_t<num_method_param, ProbClosures>, no_diffusive_t, no_source_t > ProbRHS;
//


//typedef rhs_dt<weno_t<ReconScheme::Teno5, ProbClosures>, viscous_t<visc_model_param, ProbClosures>, no_source_t > ProbRHS;



// define type of wall and EBM class

//typedef adiabatic_wall_t<ProbClosures> TypeWall;
typedef isothermal_wall_t<wall_param,ProbClosures> TypeWall; //isothermal 300 k


typedef ebm_t<TypeWall,wall_param,ProbClosures> ProbEB;



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

  // final state -- freestream
  rhot =  prob_parm.rho_oo;
  u[0] =  prob_parm.u_oo; u[1] =  prob_parm.v_oo;
  eint =  prob_parm.eint_oo;

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

  Real Rmax = 1.5;

  // refine close to square
  tagfab(i,j,k) = (Rad2 < Rmax*Rmax);
 
}
////////////////////////////////////////////////////////////////////////////////




///////////////////////////////////////////////////////////////////////////////

} // namespace PROB
#endif
