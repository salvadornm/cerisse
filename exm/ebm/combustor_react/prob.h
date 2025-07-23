// euler equations except at walls

#ifndef CNS_PROB_H_
#define CNS_PROB_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_ParmParse.H>

#include <PelePhysics.H>
#include <ReactorBase.H>

#include "Closures.h"
#include "RHS.h"

#if CNS_USE_EB    
#include <ebm.h>
#include <walltypes.h>
#endif

#include <bc_types.h>

// NTU Combustor-type for demonstration purposes
// created by S Dupre and S Navarro-Martinez (2025)


using namespace amrex;

namespace PROB {

typedef closures_dt<indicies_t, visc_suth_t,cond_suth_t ,  multispecies_pele_gas_t<indicies_t>> ProbClosures;

// problem parameters 

struct ProbParm {  
  // inflow state for initialisation
  //const Real p_inflow    = pres_atm2si; //[Pa] inflow pressure (1 atm)  
  const Real T_inflow    = 291.7;  //[K]  
  Real Y_inflow[NUM_SPECIES] = {1.0};
  ProbClosures pp_pc;

  ProbParm () {
  #if USE_PELEPHYSICS
  Y_inflow[H_ID] = 0; 
  Y_inflow[H2_ID] = 0.013; 
  Y_inflow[O_ID] = 0; 
  Y_inflow[OH_ID] = 0; 
  Y_inflow[H2O_ID] = 0; 
  Y_inflow[O2_ID] = 0.227; 
  Y_inflow[HO2_ID] = 0; 
  Y_inflow[H2O2_ID] = 0;  
  Y_inflow[N2_ID] = 0.76; 
  Y_inflow[AR_ID] = 0; 
  Y_inflow[HE_ID] = 0; 
  Y_inflow[CO_ID] = 0; 
  Y_inflow[CO2_ID] = 0; 
  #endif

  pp_pc.PYT2R(p_0, Y_inflow, T_0, rho_0);
  pp_pc.RYP2E(rho_0, Y_inflow, p_0, eint_0); 
  }

  // compute density and internal energy
  const Real Q = 8.2994; // volumetric flow rate [kg /m2 s]
  
  // inside combustor state/exit
  const Real p_0     = pres_atm2si; //[Pa] inflow pressure (1 atm) 
  const Real p_bc    = p_0;
  const Real T_0     = 285.5;  //[K]  
  Real rho_0, eint_0;

  const Real vel_0[3]= {0.0,0.0,0.0}; // array of inside velocity [m/s]

  // geometry auxiliary
  Real zin = 0.01;
  Real zexit = 0.135 ;//0.135;

};

// numerical method parameters
/** 
struct const_viscparm_t {

  public:

  //static constexpr bool dissipation = true;         // no dissipation
  //static constexpr int  order = 4;                  // order numerical scheme   
  //static constexpr Real C2skew=0.5,C4skew=0.016;   // Skew symmetric default
  static constexpr Real viscosity = 1.846e-5;
  static constexpr Real conductivity = 0.02624;
};
*/

struct viscous_param_t {

  public:
  static constexpr int order = 2;                  // order numerical scheme   

  static constexpr bool use_LES= false;

};

struct wall_param {

  public:


  static constexpr Real Twall = 285.5;                // wall temperature (if isothermal used)
  static constexpr bool solve_diffwall = true;      // solve viscous effects at walls

};


template <typename cls_t > class user_source_t;

//typedef closures_dt<indicies_t, transport_const_t<const_viscparm_t>, calorifically_perfect_gas_t<indicies_t>> ProbClosures;

// define nuemrical scheme comment/uncomment to set up 
//typedef rhs_dt<weno_t<ReconScheme::Teno5, ProbClosures>, viscous_t<viscous_param_t, ProbClosures>, user_source_t<ProbClosures> > ProbRHS;
//typedef rhs_dt<weno_t<ReconScheme::Teno5, ProbClosures>, viscous_t<viscous_param_t, ProbClosures>, no_source_t> ProbRHS;

//typedef rhs_dt<riemann_t<false, ProbClosures>, viscous_t<viscous_param_t, ProbClosures>, user_source_t<ProbClosures> > ProbRHS;
typedef rhs_dt<weno_t<ReconScheme::WenoZ5, ProbClosures>, viscous_t<viscous_param_t, ProbClosures>, reactor_source_t<user_source_t<ProbClosures>,ProbClosures >> ProbRHS;
//typedef rhs_dt<weno_t<ReconScheme::Teno5, ProbClosures>, viscous_t<viscous_param_t, ProbClosures>, reactor_source_t<user_source_t<ProbClosures>,ProbClosures >> ProbRHS;
//typedef rhs_dt<riemann_t<false, ProbClosures>, viscous_t<viscous_param_t, ProbClosures>, reactor_source_t<user_source_t<ProbClosures>,ProbClosures > > ProbRHS;



// define type of wall and EBM class
#if CNS_USE_EB    
typedef adiabatic_wall_t<ProbClosures> TypeWall;
typedef ebm_t<TypeWall,wall_param,ProbClosures> ProbEB;
#endif

typedef manual_bc_t<ProbClosures> GlobalBC;


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

  // local vars
  Real rhot,eint,u[3],y_sp[NUM_SPECIES];

  for (int n=0; n < NUM_SPECIES; n++) {
    y_sp[n] = prob_parm.Y_inflow[n];
  }
  rhot =  prob_parm.rho_0;    
  for(int idim=0;idim < AMREX_SPACEDIM;idim++) {u[idim]=prob_parm.vel_0[idim];}
  eint =  prob_parm.eint_0;

  /**
  const Real r = sqrt(x*x + y*y);
  const Real z_low = 0.048;
  const Real z_hi = 0.060;
  Real rho_spark;
  if (r <= 0.0095 && z >= z_low && z <= z_hi)
  {   
    Real Tspark = 1000.0;
    cls.PYT2R(prob_parm.p_0,y_sp, Tspark, rho_spark);
    cls.RYP2E(rho_spark, y_sp, prob_parm.p_0, eint); 
  }
  */

  Real Tspark = 1000;
  cls.PYT2R(prob_parm.p_0,y_sp, Tspark, rhot);
  cls.RYP2E(rhot, y_sp, prob_parm.p_0, eint); 
  
  
  Real kin = Real(0.5) * rhot * (u[0] * u[0] + u[1] * u[1] + u[2]*u[2]);
  //state(i, j, k, cls.URHO) = rhot;
  state(i, j, k, cls.UMX)  = rhot * u[0];
  state(i, j, k, cls.UMY)  = rhot * u[1];
  state(i, j, k, cls.UMZ)  = rhot * u[2];  
  state(i, j, k, cls.UET)  = rhot*eint + kin;    
  for (int n = 0; n < NUM_SPECIES; ++n) {
    state(i, j, k, cls.UFS + n) = rhot * y_sp[n];
  }
}
/////////////////////////////// BC /////////////////////////////////////////////
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
bcnormal(const Real x[AMREX_SPACEDIM], Real dratio, const Real s_int[ProbClosures::NCONS],
         const Real s_refl[ProbClosures::NCONS], Real s_ext[ProbClosures::NCONS], const int idir,
         const int sgn, const Real time, GeometryData const & /*geomdata*/,
         ProbClosures const &closures, ProbParm const &prob_parm) {

  const int face = (idir+1)*sgn; // +/-1 (1D) +/- 2 (2D) +/- 3 (3D)

  switch(face)
  {
    case  3:  // LEFT
	    {                  
      GlobalBC::bc_inlet_fixmassflow(0.0,0.0,1.0,&closures,
        prob_parm.Q,prob_parm.T_inflow,prob_parm.Y_inflow, s_int, s_ext);  
      break;
      }
    case  2:  // SOUTH
      GlobalBC::bc_fixP(0.0,1.0,0.0,&closures,prob_parm.p_bc, s_int, s_ext); 
      break;
    case  1:  // WEST
      GlobalBC::bc_fixP(1.0,0.0,0.0,&closures,prob_parm.p_bc, s_int, s_ext); 
      break;
    case -1:  // EAST
      GlobalBC::bc_fixP(-1.0,0.0,0.0,&closures,prob_parm.p_bc, s_int, s_ext);  
      break;
    case -2:  // NORTH
      GlobalBC::bc_fixP(-1.0,0.0,0.0,&closures,prob_parm.p_bc, s_int, s_ext); 
      break;
    case -3:   //RIGHT 
      {
      GlobalBC::bc_fixP(0.0,0.0,-1.0,&closures,prob_parm.p_bc, s_int, s_ext);  

      //GlobalBC::bc_subsonic_outflow_fixP(0.0,0.0,-1.0,&closures,prob_parm.p_0, s_int, s_ext);  
      break;  
      }
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
  Real r = sqrt(x*x + y*y);

  bool refine = false;
  
  // // refine exit of injector
  // refine= (z > 0.035) && (z < 0.07);
  // // refine close to exit  (avoid corner problem)
  // refine= (z > 0.13) || refine;


  //const int URHO= ProbClosures::URHO;

  // // compoute | d rho | normalised with rho
  // Real o_over_rhot = Real(1.0)/sdatafab(i,j,k,URHO);

  // Real drhox = Math::abs(sdatafab(i+1,j,k,URHO) - sdatafab(i-1,j,k,URHO));
  // Real drhoy = Math::abs(sdatafab(i,j+1,k,URHO) - sdatafab(i,j-1,k,URHO));
  // Real drhoz = Math::abs(sdatafab(i,j,k+1,URHO) - sdatafab(i,j,k-1,URHO));

  // Real gradrho= Real(0.5)*sqrt(drhox*drhox+drhoy*drhoy)*o_over_rhot;        


 switch (level)
  {
    case 0:
      //refine= (z > 0.035) && (z < 0.07);
      refine = (z < prob_parm.zexit); 
      //refine = (z < prob_parm.zexit) && (r < 0.025) ;    // refine combustor    
      break;
    case 1:
      refine= (z > 0.035) && (z < 0.07);    
      // refine = (z < prob_parm.zexit); 
      break;
    case 2:
      // refine= (z > 0.035) && (z < 0.07);    
      break;  
      
    default:

    break;
   }

  tagfab(i,j,k) = refine;


}
///////////////////////////////SOURCE TERM /////////////////////////////////////
template <typename cls_t>
class user_source_t {
  public:


  // to use as a user source term, the function name must be rsrc:
  // void inline src(const Geometry& geomdata, const amrex::MFIter &mfi,
  //                 const amrex::Array4<const amrex::Real> &prims,
  //                 const amrex::Array4<amrex::Real> &rhs, const cls_t *cls_d,
  //                 amrex::Real dt){

 
  // to use as a user source term after calling reactor, the function name must be rsrc:
  void inline rsrc(const Geometry& geomdata, const amrex::MFIter &mfi,
                  const amrex::Array4<const amrex::Real> &prims,
                  const amrex::Array4<amrex::Real> &rhs, const cls_t *cls_d,
                  amrex::Real dt){

    const Box& bxg = mfi.tilebox();
    // const Box& bxg = mfi.growntilebox(cls_t::NGHOST);
    const Real *prob_lo = geomdata.ProbLo();
    const Real *dx = geomdata.CellSize();

    ProbParm const prob_parm;
    const auto& cls = *cls_d;


    amrex::ParallelFor(bxg,
      [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
      
      const Real x = prob_lo[0] + (i + Real(0.5)) * dx[0];
      const Real y = prob_lo[1] + (j + Real(0.5)) * dx[1];
      const Real z = prob_lo[2] + (k + Real(0.5)) * dx[2];
      const Real r = sqrt(x*x + y*y + (z-prob_parm.zexit)*(z-prob_parm.zexit));
      
      // pressure relax  if P > P0  P drops and to keep T constant rho drops
      const Real tau_relax = 50.0*dt; 
      const Real coef =dt/tau_relax;

      const Real dP = (prob_parm.p_bc- prims(i,j,k,cls.QPRES))*coef;
    
      const Real T = prims(i,j,k,cls.QT); // dT =0
      const Real rho  = prims(i, j, k, cls.QRHO);
     
      Real drho = 0.0; Real Y[NUM_SPECIES] = {1.0};
      #if USE_PELEPHYSICS
      for (int sp = 0; sp < NUM_SPECIES; sp ++) {
        Y[sp] = prims(i, j, k, cls.QFS + sp);
      }
      #endif 

      cls.PYT2R(dP,Y,T,drho);  
      Real drhodt = drho/dt;
      //Real drhodt = drho;
      Real kin = 0.5*(prims(i,j,k,cls.QU)*prims(i,j,k,cls.QU) + prims(i,j,k,cls.QV)*prims(i,j,k,cls.QV)
                    + prims(i,j,k,cls.QW)*prims(i,j,k,cls.QW));
      Real Et  = prims(i,j,k,cls.QEINT) + kin;
    
      bool buffer = (z > prob_parm.zexit) && (r > 0.035);

      if (buffer){        
        rhs(i,j,k,cls.UMX) += prims(i,j,k,cls.QU)*drhodt;
        rhs(i,j,k,cls.UMY) += prims(i,j,k,cls.QV)*drhodt; 
        rhs(i,j,k,cls.UMZ) += prims(i,j,k,cls.QW)*drhodt;
        rhs(i,j,k,cls.UET) += Et*drhodt;
        for (int sp = 0; sp < NUM_SPECIES; sp++) {
          rhs(i,j,k, cls.UFS + sp) += prims(i,j,k, cls.QFS + sp) * drhodt;  
        }                        
      }


      });

  };
};

} // namespace PROB

#endif
