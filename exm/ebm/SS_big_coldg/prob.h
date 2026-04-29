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
#include <numbers>
#include <cmath>


// NTU Combustor-type for demonstration purposes
// created by S Dupre and S Navarro-Martinez (2025)


struct NTNUglobalnumbers {  
  static constexpr Real pressure = pres_atm2si; // [Pa] inflow pressure (1 atm)  
  static constexpr Real Temperature_inflow = 298.0; //[K] 
  // geometry 
  static constexpr Real zinlet = 0.045;  // inlet combustor
  static constexpr Real zexit  = 0.135 ; // outlet combustor  
  static constexpr Real radius = 0.022 ; // Radius combustor    
};

/*
struct TFMparm {
  amrex::Real F0      = 1.0;
  amrex::Real SLO     = 2.0;
  amrex::Real deltaLO = 0.0001;
  amrex::Real Ret     = 10.0;
  int sensor_id       = 0;     // set later from idx_t::QT
  amrex::Real c1      = 2200.0;
  amrex::Real c0      = 298.0;
  amrex::Real beta    = 1.0;
};
*/

using namespace amrex;

namespace PROB {


// LES closures

struct LESparm {
  // Smagorinsky constant
  static constexpr Real Cs = 0.1; // Smag original 0.17 (off)
  static constexpr int order = 2; // order of the numerical scheme for LES
  static constexpr Real Scsgs = 0.4 ; //0.4 // turbulent Schmidt number
  static constexpr Real Pr_o_Prsgs = 1.0; //0.1 // turbulent Prandtl number  
  static constexpr bool fixDelta = false; // use fixed filter width
};


// Closure Index+Thermodynamics + Transport + LES + ATF

typedef closures_dt<indicies_stat_t, transport_Pele_t , multispecies_pele_gas_t<indicies_t>,
                    Smagorinsky_t<LESparm,indicies_t>, TFM_t<LESparm,indicies_t>>ProbClosures;
//typedef closures_dt<indicies_t, transport_Pele_t , multispecies_pele_gas_t<indicies_t>> ProbClosures;

// problem parameters 

struct ProbParm {  
  // inflow state for initialisation
  //const Real p_inflow    = pres_atm2si; //[Pa] inflow pressure (1 atm)  
  const Real T_inflow    = 298.0; //[K]  
  Real Y_inflow[NUM_SPECIES] = {0.0};
  Real Y_burn[NUM_SPECIES] = {0.0};
  ProbClosures pp_pc;
  NTNUglobalnumbers combustor;

  const Real Tburn = 1426.8; //[K] burned gas temperature


  // inflow
  const Real p_0     = combustor.pressure; //[Pa] inflow pressure (1 atm) 
  const Real T_0     = combustor.Temperature_inflow;  //[K]  

  // compute density and internal energy
  const Real Q = 8.665; //  flow m3/s

  // aux variables for initialisation
  Real rho_0, eint_0;


  ProbParm () {
#if USE_PELEPHYSICS
  Y_inflow[H_ID]    = 0; 
  Y_inflow[H2_ID]   = 0.011608; 
  Y_inflow[O_ID]    = 0; 
  Y_inflow[OH_ID]   = 0; 
  Y_inflow[H2O_ID]  = 0; 
  Y_inflow[O2_ID]   = 0.2303; 
  Y_inflow[HO2_ID]  = 0; 
  Y_inflow[H2O2_ID] = 0;  
  Y_inflow[N2_ID]   = 0.7581; 
   
  // equivalence ratio = 0.4
  Y_burn[H_ID]    = 1.4351e-10; //4.0955e-08;//9.125e-06; 
  Y_burn[H2_ID]   = 2.2792e-08; //1.742e-06;//7.2983e-05; 
  Y_burn[O_ID]    = 2.8586e-07; //1.2649e-05;//0.00068717; 
  Y_burn[OH_ID]   = 3.2613e-05; //0.00042606;//0.0047711; 
  Y_burn[H2O_ID]  = 0.10371;    //0.14176;//0.12602; 
  Y_burn[O2_ID]   = 0.13816;    // 0.10298;//0.11253; 
  Y_burn[HO2_ID]  = 1.1152e-07; //8.168e-07;//5.6389e-06; 
  Y_burn[H2O2_ID] = 8.6331e-09; //5.5039e-08;//3.3351e-07;  
  Y_burn[N2_ID]   = 0.7581;     //0.75481;//0.7559 ; 
  
#endif

  pp_pc.PYT2R(p_0, Y_inflow, T_0, rho_0);
  pp_pc.RYP2E(rho_0, Y_inflow, p_0, eint_0); 
  }
  
  // 
  const Real T_b  = Tburn;
  const Real T_u  = T_0;
  
  const Real vel_0[3]= {0.0,0.0,0.0}; // array of inside velocity [m/s]

};

struct IgnitionParm {
  static constexpr bool ignite = true;
  // Time Window (in seconds)
  static constexpr int dt_min = 0;
  static constexpr int dt_max = 10000;
  // Spatial Parameters
  static constexpr Real x0 = 0;
  static constexpr Real y0 = 0;
  static constexpr Real z0 = 0.05;
  static constexpr Real x_sd = 0.00175;
  static constexpr Real y_sd = 0.00175;
  static constexpr Real z_sd = 0.00175;
  static constexpr Real nsigma = 3.0;
  static constexpr Real x_min = x0 - nsigma * x_sd;
  static constexpr Real x_max = x0 + nsigma * x_sd;
  static constexpr Real y_min = y0 - nsigma * y_sd;
  static constexpr Real y_max = y0 + nsigma * y_sd;
  static constexpr Real z_min = z0 - nsigma * z_sd;
  static constexpr Real z_max = z0 + nsigma * z_sd;
  // Tanh scaling parameters
  static constexpr Real p = 0.99;
  static constexpr Real k = 2.646652412362245 / (dt_max - dt_min); // atanh(0.99) ~ 2.64
  // Spark parameters
  static constexpr Real Tmax = 1800;
};

// numerical method parameters 
struct skewparm_t {

  public:

  static constexpr bool dissipation = true;         // no dissipation
  static constexpr int  order = 4;                  // order numerical scheme   (2 or 4)
  static constexpr Real C2skew=0.1,C4skew=0.016;    // Skew symmetric default  (0.1 0.016)
};


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


template <typename cls_t, typename cls_ignite > class user_source_t;

// USED
//typedef rhs_dt<riemann_t<false, ProbClosures>, viscousLES_t<user_source_t<ProbClosures>, ProbClosures>, reactor_sourceLES_t<user_source_t<ProbClosures>,ProbClosures >> ProbRHS;
typedef rhs_dt<skew_t<skewparm_t, ProbClosures>, viscousLES_t<user_source_t<ProbClosures,IgnitionParm>, ProbClosures>, reactor_sourceLES_t<user_source_t<ProbClosures,IgnitionParm>,ProbClosures >> ProbRHS;
//typedef rhs_dt<weno_t<ReconScheme::Teno5, ProbClosures>, viscousLES_t<user_source_t<ProbClosures>, ProbClosures>, reactor_sourceLES_t<user_source_t<ProbClosures>,ProbClosures >> ProbRHS;


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

  //Real x = prob_lo[0] + (i + Real(0.5)) * dx[0];
  //Real y = prob_lo[1] + (j + Real(0.5)) * dx[1];
  //Real z = prob_lo[2] + (k + Real(0.5)) * dx[2];

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

  // put burn conditions inside combustor
  bool comb_init = false;
  if (comb_init)
  {
    // if (z > 0.035)
    // { 
      for (int n=0; n < NUM_SPECIES; n++) {
        y_sp[n] = prob_parm.Y_burn[n];
      }     
      cls.PYT2R(prob_parm.p_0,y_sp, prob_parm.Tburn, rhot);
      cls.RYP2E(rhot, y_sp, prob_parm.p_0, eint);
    // }
    // else
    // {
    //   u[2] = prob_parm.Q/prob_parm.rho_0;
    // }


  }

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

  //snm comment for 0-gradient BCs (to be used for debugging)
  //for (int n=0; n < ProbClosures::NCONS; n++) {
  //  s_ext[n] = s_int[n];
 // }
 // return;


  switch(face)
  {
    case  3:  // LEFT/BOTTOM  z
	    {                  
      GlobalBC::bc_inlet_fixmassflow(0.0,0.0,1.0,&closures,
        prob_parm.Q,prob_parm.T_inflow,prob_parm.Y_inflow, s_int, s_ext);  
        
      break;
      }
    case  2:  // SOUTH        y  
      GlobalBC::bc_fixP(0.0,1.0,0.0,&closures,prob_parm.p_0, s_int, s_ext); 
      break;
    case  1:  // WEST         x  
      GlobalBC::bc_fixP(1.0,0.0,0.0,&closures,prob_parm.p_0, s_int, s_ext); 
      break;
    case -1:  // EAST         x  
      GlobalBC::bc_fixP(-1.0,0.0,0.0,&closures,prob_parm.p_0, s_int, s_ext);  
      break;
    case -2:  // NORTH        y 
      GlobalBC::bc_fixP(0.0,-1.0,0.0,&closures,prob_parm.p_0, s_int, s_ext); 
      break;
    case -3:   // RIGHT/TOP   z
      {
      GlobalBC::bc_fixP(0.0,0.0,-1.0,&closures,prob_parm.p_0, s_int, s_ext);  

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

  NTNUglobalnumbers const combustor;

  Real x = prob_lo[0] + (i + Real(0.5)) * dx[0];
  Real y = prob_lo[1] + (j + Real(0.5)) * dx[1];
  Real z = prob_lo[2] + (k + Real(0.5)) * dx[2];
  Real r = sqrt(x*x + y*y);

  bool refine = false;
  

  //  U->Q 
  Real Q[ProbClosures::NPRIM],U[ProbClosures::NCONS];
  for (int n = 0; n < ProbClosures::NCONS; ++n) {
    U[n] = sdatafab(i,j,k,n);
  }
  auto thermo = ProbClosures::multispecies_pele_gas_t();
  thermo.cons2prims_point(U,Q);
 

  // use H2O as refinement
  //const int UREF= ProbClosures::UFS + H2O_ID;  
  // Real drhox =( sdatafab(i+1,j,k,UREF) - sdatafab(i-1,j,k,UREF) );
  // Real drhoy =( sdatafab(i,j+1,k,UREF) - sdatafab(i,j-1,k,UREF) );
  // Real drhoz =( sdatafab(i,j,k+1,UREF) - sdatafab(i,j,k-1,UREF) );
  // Real gradrho= Real(0.5)*sqrt(drhox*drhox+drhoy*drhoy);
  // refine = grad_rho > 0.1       
  
  // progress variable
  //Real c =  (Q[ProbClosures::QT]-prob_parm.T_u)/(prob_parm.T_b - prob_parm.T_u);


 switch (level)
  {
    case 0:      
      //refine = (z < combustor.zexit) && (r < 1.1*combustor.radius) ;   
      refine = (z >= combustor.zexit) && (z < 0.14) ||  (z < combustor.zexit) && (r < 1.2*combustor.radius);
      break;
    case 1:
      refine = (z < combustor.zexit) && (r < 1.1*combustor.radius);   
      //refine = (z < 0.08) && (r < 0.025);  
      break;
    case 2:
      refine = (z > 0.03) && (z < 0.06) && (r < 0.015);
      break;  
    case 3:  
      //refine =  (c < 0.7) && (c > 0.3);
      break;      
    default:

    break;
  }

  tagfab(i,j,k) = refine;


}
///////////////////////////////SOURCE TERM /////////////////////////////////////
template <typename cls_t, typename cls_ignite>
class user_source_t {
  public:

  // ATF options
  bool static constexpr use_ATF = true; // use adaptive thickening factor
  static constexpr int ATF_model = 1;   // 1: Classic  Colin/Charlette , 2: Rathore transformation
  //

  //...
  
  // viscous LES options
  static constexpr bool use_LES=true;           // use LES model in viscous term
  static constexpr int  order = 2;               // order numerical scheme   

  // compute chemistry
  bool static constexpr do_reactions = true;  // <<<<<<<<<<<<<<<<<<<<<<<<<
  bool static constexpr mask_cells_boundary = true; // avoid compute reactions in cells partially covered
  //

  inline static int src_dt = 0;

  static void inline increase_src_dt(){
    src_dt++;
  }

  // to use as a user source term, the function name must be src:
  // void inline src(const Geometry& geomdata, const amrex::MFIter &mfi,
  //                 const amrex::Array4<const amrex::Real> &prims,
  //                 const amrex::Array4<amrex::Real> &rhs, const cls_t *cls_d,
  //                 amrex::Real dt){

 
  // to use as a user source term after calling reactor, the function name must be rsrc:
  void inline rsrc(const Geometry& geomdata, const amrex::MFIter &mfi,
                  const amrex::Array4<const amrex::Real> &prims,
                  const amrex::Array4<amrex::Real> &rhs, const cls_t *cls_d,
                  amrex::Real dt, amrex::Real real_time, const int inc_dt){

    const Box& bxg = mfi.tilebox();
    
    // use device-friendly arrays instead of pointers
    //  const Real *prob_lo = geomdata.ProbLo();
    //  const Real *dx = geomdata.CellSize();
    auto prob_lo = geomdata.ProbLoArray();
    auto dx     = geomdata.CellSizeArray();
                
    NTNUglobalnumbers const combustor;

    const Real tau_relax = 5.e-6;
    const Real coef =dt/tau_relax;
 
    const Real zbuffer  = 0.16;
    const Real p_0      = combustor.pressure;

    // ------------------------------------------------------------------------------------
    amrex::ParallelFor(bxg,
      [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
      
      // coordinates 
      const Real x = prob_lo[0] + (i + Real(0.5)) * dx[0];
      const Real y = prob_lo[1] + (j + Real(0.5)) * dx[1];
      const Real z = prob_lo[2] + (k + Real(0.5)) * dx[2];
      const Real r = sqrt(x*x + y*y + (z-combustor.zexit)*(z-combustor.zexit));

      // get values from primitives
      const Real pres = prims(i,j,k,cls_t::QPRES);
      const Real T    = prims(i,j,k,cls_t::QT   ); 
      const Real rho  = prims(i,j,k,cls_t::QRHO );     
      Real Y[NUM_SPECIES]; for (int sp = 0; sp < NUM_SPECIES; sp ++) { Y[sp] = prims(i, j, k, cls_t::QFS + sp);}
      // u v w      
      const Real u  = prims(i,j,k,cls_t::QU );
      const Real v  = prims(i,j,k,cls_t::QV );
      const Real w  = prims(i,j,k,cls_t::QW );
      // Kinetic and Total energy (per density unit)  [J/rho] 
      const Real Et  = prims(i,j,k,cls_t::QEINT) + 0.5*(u*u + v*v + w*w);

      // pressure relaxation towards P0   
      const Real dP = (p_0- pres)*coef; Real drho;
      cls_d->PYT2R(dP,Y,T,drho);      // associated density change
      const Real drhodt = drho/dt;  // rate of change
   
      // buffer zone where relaxation applies 
      const bool buffer = (z > combustor.zexit) && (r > 0.06);

      if (buffer){        
 
        rhs(i,j,k,cls_t::UMX) += u*drhodt;
        rhs(i,j,k,cls_t::UMY) += v*drhodt; 
        rhs(i,j,k,cls_t::UMZ) += w*drhodt;
        rhs(i,j,k,cls_t::UET) += Et*drhodt;
        for (int sp = 0; sp < NUM_SPECIES; sp++)  rhs(i,j,k, cls_t::UFS + sp) += Y[sp] * drhodt;  
      }
      
      if constexpr (!cls_ignite::ignite) return;

      if (inc_dt > cls_ignite::dt_min && inc_dt < cls_ignite::dt_max){

        if  (x >= cls_ignite::x_min && x <= cls_ignite::x_max && y >= cls_ignite::y_min && y <= cls_ignite::y_max && z >= cls_ignite::z_min && z <= cls_ignite::z_max)
        {
          const Real i_mean = Real(std::round(float((cls_ignite::x0 - prob_lo[0]) / dx[0] - Real(0.5)))); 
          const Real j_mean = Real(std::round(float((cls_ignite::y0 - prob_lo[1]) / dx[1] - Real(0.5)))); 
          const Real k_mean = Real(std::round(float((cls_ignite::z0 - prob_lo[2]) / dx[2] - Real(0.5)))); 
          const Real i_sd = Real(std::round(std::abs(float(cls_ignite::x_sd / dx[0]))));
          const Real j_sd = Real(std::round(std::abs(float(cls_ignite::y_sd / dx[1]))));
          const Real k_sd = Real(std::round(std::abs(float(cls_ignite::z_sd / dx[2]))));
        
          Real rho_ignite, eint_ignite;
          cls_d->PYT2R(pres, Y, cls_ignite::Tmax, rho_ignite);
          cls_d->RYP2E(rho_ignite, Y, pres, eint_ignite);
          eint_ignite = std::max((Real)0, eint_ignite - prims(i,j,k, cls_t::QEINT));
          eint_ignite /= dt;

          // spatial distribution of nrj source
          Real gaussian_space = std::exp( -0.5 * ( (i - i_mean)*(i - i_mean)/(i_sd*i_sd) + (j - j_mean)*(j - j_mean)/(j_sd*j_sd) + 
                                    (k - k_mean)*(k - k_mean)/(k_sd*k_sd)));

          // temporal distribution of nrj source
          Real tanh_time = Real(std::tanh(float(cls_ignite::k * inc_dt)));

          // total nrj source 
          Real tanh_total = rho * gaussian_space * tanh_time * eint_ignite;
  
          rhs(i,j,k,cls_t::UET) += tanh_total;
        }
      }

    });

  };
};

} // namespace PROB

#endif
