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


using namespace amrex;

namespace PROB {


// LES closures

struct LESparm {
  // Smagorinsky constant
  static constexpr Real Cs = 0.0; // Smag original 0.17 (off)
  static constexpr int order = 2; // order of the numerical scheme for LES
  static constexpr Real Scsgs = 0.4 ; //0.4 // turbulent Schmidt number
  static constexpr Real Pr_o_Prsgs = 1.0; //0.1 // turbulent Prandtl number  
  static constexpr bool fixDelta = false; // use fixed filter width
};


typedef closures_dt<indicies_stat_t, transport_Pele_t , multispecies_pele_gas_t<indicies_t>, Smagorinsky_t<LESparm,indicies_t>>ProbClosures;
//typedef closures_dt<indicies_t, transport_Pele_t , multispecies_pele_gas_t<indicies_t>> ProbClosures;

// problem parameters 

struct ProbParm {  
  // inflow state for initialisation
  //const Real p_inflow    = pres_atm2si; //[Pa] inflow pressure (1 atm)  
  const Real T_inflow    = 298.0; //[K]  
  Real Y_inflow[NUM_SPECIES] = {0.0};
  Real Y_burn[NUM_SPECIES] = {0.0};
  ProbClosures pp_pc;

  const Real Tburn = 2207.0; //[K] burned gas temperature

  ProbParm () {
  #if USE_PELEPHYSICS
  Y_inflow[H_ID] = 0; 
  Y_inflow[H2_ID] = 0.014468; 
  Y_inflow[O_ID] = 0; 
  Y_inflow[OH_ID] = 0; 
  Y_inflow[H2O_ID] = 0; 
  Y_inflow[O2_ID] = 0.22963; 
  Y_inflow[HO2_ID] = 0; 
  Y_inflow[H2O2_ID] = 0;  
  Y_inflow[N2_ID] = 0.7559; 
   
  //
  Y_burn[H_ID] = 9.125e-06; 
  Y_burn[H2_ID] = 7.2983e-05; 
  Y_burn[O_ID] = 0.00068717; 
  Y_burn[OH_ID] = 0.0047711; 
  Y_burn[H2O_ID] = 0.12602; 
  Y_burn[O2_ID] = 0.11253; 
  Y_burn[HO2_ID] = 5.6389e-06; 
  Y_burn[H2O2_ID] = 3.3351e-07;  
  Y_burn[N2_ID] = 0.7559 ; 
  
  #endif

  pp_pc.PYT2R(p_0, Y_inflow, T_0, rho_0);
  pp_pc.RYP2E(rho_0, Y_inflow, p_0, eint_0); 
  }

  // compute density and internal energy
  const Real Q = 8.665; //  flow m3/s
  
  // inflow
  const Real p_0     = pres_atm2si; //[Pa] inflow pressure (1 atm) 
  const Real T_0     = 298.0;  //[K]  
  
  // 
  const Real T_b  = Tburn;
  const Real T_u  = T_0;
  

  Real rho_0, eint_0;

  const Real vel_0[3]= {0.0,0.0,0.0}; // array of inside velocity [m/s]

  // geometry auxiliary
  Real zin = 0.01;
  Real zexit = 0.135 ;//0.135;

};

// spark parametrs
struct SparkParm{  
  const Real t0 = 0.0085;          // spark time 
  const Real x0 = 0.0;          // spark position
  const Real y0 = 0.0;
  const Real z0 = 0.08;     
  const Real a  = 4.0*std::sqrt(std::log(10));
  const Real Tmax    = 3000.0;
  const Real energy  = 2000.0e-3; // 100 mJ
  const Real pi = 3.14159265359;
  //const Real ds   = std::sqrt(a/pi)*(energy);
  const Real ds    = 4.0e-3;    // 5  mm
  const Real dt    = 1.0e-3;    // 0.5 ms
  const Real dt2   = dt*dt;     // 
  const Real ds2   = ds*ds;     // 
  const Real o_Volt = 0.25/(pi*pi*ds*ds*ds*dt);
  const Real Cp0    =  1224; // [J/(kg K) assumed room temperature and phi=0.5
  const Real  T0    = 300;
  const Real  rho0  = 0.97 ;
  const Real  ds0   = sqrt(a/pi)*std::pow(energy/(rho0*Cp0*(Tmax-T0)), 1.0/3.0);
  // ds size that corresponds to a maximum temperature of Tmax, with given energy
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


template <typename cls_t > class user_source_t;

// USED
//typedef rhs_dt<weno_t<ReconScheme::WenoZ5, ProbClosures>, viscousLES_t<user_source_t<ProbClosures>, ProbClosures>, reactor_sourceLES_t<user_source_t<ProbClosures>,ProbClosures >> ProbRHS;
//typedef rhs_dt<skew_t<skewparm_t, ProbClosures>, viscousLES_t<user_source_t<ProbClosures>, ProbClosures>, reactor_sourceLES_t<user_source_t<ProbClosures>,ProbClosures >> ProbRHS;

//typedef rhs_dt<skew_t<skewparm_t, ProbClosures>, no_diffusive_t, no_source_t > ProbRHS;

typedef rhs_dt<skew_t<skewparm_t, ProbClosures>, viscousLES_t<user_source_t<ProbClosures>, ProbClosures>, no_source_t > ProbRHS;


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

  // put burn conditions inside combustor
  bool comb_init = true;
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

  // snm comment for 0-gradient BCs (to be used for debugging)
  // for (int n=0; n < ProbClosures::NCONS; n++) {
  //   s_ext[n] = s_int[n];
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
  Real x = prob_lo[0] + (i + Real(0.5)) * dx[0];
  Real y = prob_lo[1] + (j + Real(0.5)) * dx[1];
  Real z = prob_lo[2] + (k + Real(0.5)) * dx[2];
  Real r = sqrt(x*x + y*y);

  bool refine = false;
  

  //  U->Q 
  // Real Q[ProbClosures::NPRIM],U[ProbClosures::NCONS];
  // for (int n = 0; n < ProbClosures::NCONS; ++n) {
  //   U[n] = sdatafab(i,j,k,n);
  // }
  // auto thermo = ProbClosures::multispecies_pele_gas_t();
  // thermo.cons2prims_point(U,Q);
 

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
      refine =(z < 0.08) && (r < 0.1);    
      break;
    case 1:
      refine = (z > 0.03) && (z < 0.08) && (r < 0.015);  
      break;
    case 2:
      refine = (z > 0.0042) && (z < 0.06) && (r < 0.01);
      break;  
    case 3:  
      //refine =  (c < 0.7) && (c > 0.3);
      break;      
    default:

    break;
  }

 // refine = true; // temp

  tagfab(i,j,k) = refine;


}
///////////////////////////////SOURCE TERM /////////////////////////////////////
template <typename cls_t>
class user_source_t {
  public:

  // ATF options
  bool static constexpr ATF = true; // use adaptive thickening factor
  static constexpr Real thickfactor = 1.0; // thickening factor
  // compute chemistry
  bool static constexpr do_reactions = true;  // <<<<<<<<<<<<<<<<<<<<<<<<<
  bool static constexpr mask_cells_boundary = true; // avoid compute reactions in cells partially covered
  //

  // viscous LES options
  static constexpr int order = 2;                  // order numerical scheme   
  static constexpr bool use_LES= false;           // use LES model in viscous term


  // to use as a user source term, the function name must be src:
  // void inline src(const Geometry& geomdata, const amrex::MFIter &mfi,
  //                 const amrex::Array4<const amrex::Real> &prims,
  //                 const amrex::Array4<amrex::Real> &rhs, const cls_t *cls_d,
  //                 amrex::Real dt){

 
  // to use as a user source term after calling reactor, the function name must be rsrc:
  void inline rsrc(const Geometry& geomdata, const amrex::MFIter &mfi,
                  const amrex::Array4<const amrex::Real> &prims,
                  const amrex::Array4<amrex::Real> &rhs, const cls_t *cls_d,
                  amrex::Real dt, amrex::Real real_time){

    const Box& bxg = mfi.tilebox();
    // const Box& bxg = mfi.growntilebox(cls_t::NGHOST);
    const Real *prob_lo = geomdata.ProbLo();
    const Real *dx = geomdata.CellSize();

    ProbParm const prob_parm;
    const auto& cls = *cls_d;

    // for spark
    //SparkParm const spark;
    //Real gauss_funt = std::exp( -0.5 * (real_time - spark.t0)*(real_time - spark.t0)/spark.dt2 );


    const Real tau_relax = 5.e-6;
    const Real coef =dt/tau_relax;
 

    // doing in this way avoids passing  whole prob_parm and only pass scalars you use 
    const Real zexit = prob_parm.zexit;
    const Real p_0    = prob_parm.p_0;

    // ------------------------------------------------------------------------------------
    amrex::ParallelFor(bxg,
      [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
      
      const Real x = prob_lo[0] + (i + Real(0.5)) * dx[0];
      const Real y = prob_lo[1] + (j + Real(0.5)) * dx[1];
      const Real z = prob_lo[2] + (k + Real(0.5)) * dx[2];
      const Real r = sqrt(x*x + y*y + (z-zexit)*(z-zexit));
      
      // pressure relax  if P > P0  away from chamber
      Real pres = prims(i,j,k,cls.QPRES);

      const Real dP = (p_0- pres)*coef;
    
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
    
      bool buffer = (z > zexit) && (r > 0.06);

      if (buffer){        
        rhs(i,j,k,cls.UMX) += prims(i,j,k,cls.QU)*drhodt;
        rhs(i,j,k,cls.UMY) += prims(i,j,k,cls.QV)*drhodt; 
        rhs(i,j,k,cls.UMZ) += prims(i,j,k,cls.QW)*drhodt;
        rhs(i,j,k,cls.UET) += Et*drhodt;
        for (int sp = 0; sp < NUM_SPECIES; sp++) {
          rhs(i,j,k, cls.UFS + sp) += prims(i,j,k, cls.QFS + sp) * drhodt;  
        }                        
      }

      // spark ignition --------------------
      // const Real rs2 = (x- spark.x0)*(x- spark.x0) + (y- spark.y0)*(y- spark.y0) + (z- spark.z0)*(z- spark.z0);
      // Real gauss_funr = std::exp( -0.5 * rs2/spark.ds2 );        
      // rhs(i,j,k,cls.UET) += spark.energy*gauss_funt*gauss_funr*spark.o_Volt ;

      
      // damps energy if T > 3000 K
      // if (T > 3000.0)
      // {
      //   const Real Ts = 3000.0;       
      //   Real Etarget = 0.0;
      //   Real rhos = 0.0;
      //   cls.PYT2R(pres,Y,Ts,rhos);  
      //   cls.RYP2E(rhos,Y,pres,Etarget);
      //   Real tau = 5.0e-6; //5e-3 slow relaxation
      //   rhs(i,j,k,cls.UET) += (rhos*Etarget - rho*prims(i,j,k,cls.QEINT))/tau;
      // }

    });

  };
};

} // namespace PROB

#endif
