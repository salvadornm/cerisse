// euler equations except at walls

#ifndef CNS_PROB_H_
#define CNS_PROB_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_ParmParse.H>

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

struct CavityParm {
  //free stream conditions at p = 1 atm
  static constexpr Real rho_inf = 1.177; //[kg/m^3]
  static constexpr Real T_inf = 300; //[K]
  static constexpr Real viscosity = 1.845e-5; //[Pa.s]
  static constexpr Real conductivity = 0.0262; //[W/K.m]
  static constexpr Real R = 287.05; //[J/kg.K]
  static constexpr Real gamma = 1.4; 
  static constexpr Real Ma = 0.6;
  static constexpr Real Re = 56.8;
  static constexpr Real theta = 4.274e-6; //[m]
  static constexpr Real x_virtual = 0.000190;
};

// LES closures
struct LESparm {
  // Smagorinsky constant
  static constexpr Real Cs = 0.0; //0.1 prev. 0.0 = 0ff
  static constexpr int order = 2; // order of the numerical scheme for LES
  static constexpr Real Scsgs = 0.4 ; //0.4 // turbulent Schmidt number
  static constexpr Real Pr_o_Prsgs = 0.9; //0.1 // turbulent Prandtl number  
  static constexpr bool fixDelta = false; // use fixed filter width
};


typedef closures_dt<indicies_stat_t, transport_const_t<CavityParm>, calorifically_perfect_gas_t<indicies_t>, Smagorinsky_t<LESparm,indicies_t>>ProbClosures;
//typedef closures_dt<indicies_t, transport_Pele_t , multispecies_pele_gas_t<indicies_t>> ProbClosures;

// problem parameters 

struct ProbParm {  
  CavityParm cavity;
  ProbClosures pp_pc;

  // inflow state for initialisation
  const Real T_0    = cavity.T_inf;
  const Real p_0     = pres_atm2si; //[Pa] inflow pressure (1 atm)   
  Real rho_0, eint_0, c_0;
  Real vel_inlet[3] = {0.0,0.0,0.0};
  Real Q;

  //geom
  Real zexit = 0;

  ProbParm () {
  pp_pc.PYT2R(p_0, nullptr, T_0, rho_0);
  pp_pc.RYP2E(rho_0, nullptr, p_0, eint_0);
  pp_pc.PYT2Cs(p_0, nullptr, T_0, c_0);
  vel_inlet[0] = cavity.Ma * c_0;
  Q = cavity.Re * cavity.viscosity / cavity.theta;
  }
};

// numerical method parameters 
struct skewparm_t {

  public:

  static constexpr bool dissipation = false;         // no dissipation
  static constexpr int  order = 6;                  // order numerical scheme   (2 or 4)
  static constexpr Real C2skew=0.1,C4skew=0.016;    // Skew symmetric default  (0.5)
};


struct viscous_param_t {

  public:
  static constexpr int order = 2;                  // order numerical scheme   
  static constexpr bool use_LES= false;
};

struct wall_param {

  public:

  static constexpr Real Twall = CavityParm::T_inf;                // wall temperature (if isothermal used)
  static constexpr bool solve_diffwall = true;      // solve viscous effects at walls

};


template <typename cls_t > class user_source_t;

// USED
//typedef rhs_dt<weno_t<ReconScheme::WenoZ5, ProbClosures>, viscousLES_t<user_source_t<ProbClosures>, ProbClosures>, reactor_sourceLES_t<user_source_t<ProbClosures>,ProbClosures >> ProbRHS;
//typedef rhs_dt<riemann_t<false, ProbClosures>, viscousLES_t<user_source_t<ProbClosures>, ProbClosures>, reactor_sourceLES_t<user_source_t<ProbClosures>,ProbClosures >> ProbRHS;
//typedef rhs_dt<skew_t<skewparm_t, ProbClosures>, viscousLES_t<user_source_t<ProbClosures>, ProbClosures>, reactor_sourceLES_t<user_source_t<ProbClosures>,ProbClosures >> ProbRHS;

// Skew Symmetric. No LES
typedef rhs_dt<skew_t<skewparm_t, ProbClosures>, viscousLES_t<user_source_t<ProbClosures>, ProbClosures>, user_source_t<ProbClosures>> ProbRHS;


// define type of wall and EBM class
#if CNS_USE_EB    
//typedef adiabatic_wall_t<ProbClosures> TypeWall;
typedef isothermal_wall_t<wall_param,ProbClosures> TypeWall;
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
  CavityParm cavity;

  Real x = prob_lo[0] + (i + Real(0.5)) * dx[0];
  Real y = prob_lo[1] + (j + Real(0.5)) * dx[1];
  Real z = prob_lo[2] + (k + Real(0.5)) * dx[2];

  // local vars
  Real rhot,eint,u[3];

  rhot =  prob_parm.rho_0;

  //free stream velocities    
  for(int idim=0;idim < AMREX_SPACEDIM;idim++) {u[idim]=prob_parm.vel_inlet[idim];}
  //Paulhausen profile
  Real x_0 = cavity.x_virtual;
  Real delta = 5.29 * sqrt((cavity.viscosity * (x + x_0)) / (cavity.rho_inf * prob_parm.vel_inlet[0]));
  if (y < delta) {u[0] = (2 * (y/delta) - 2 * (y/delta)*(y/delta)*(y/delta) + (y/delta)*(y/delta)*(y/delta)*(y/delta)) * prob_parm.vel_inlet[0];}
  if (y <= 0) {u[0] = 0;}

  eint =  prob_parm.eint_0;

  Real kin = Real(0.5) * rhot * (u[0] * u[0] + u[1] * u[1] + u[2]*u[2]);
  state(i, j, k, cls.URHO) = rhot;
  state(i, j, k, cls.UMX)  = rhot * u[0];
  state(i, j, k, cls.UMY)  = rhot * u[1];
  state(i, j, k, cls.UMZ)  = rhot * u[2];  
  state(i, j, k, cls.UET)  = rhot*eint + kin;    
  
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
    case  1:  // west x
	    {                  
      GlobalBC::bc_inlet_fixmassflow(1.0,0.0,0.0,&closures,
        prob_parm.Q,prob_parm.T_0, nullptr, s_int, s_ext);  
      break;
      }
    case -1:  // EAST         x  
      GlobalBC::bc_fixP(-1.0,0.0,0.0,&closures,prob_parm.p_0, s_int, s_ext);  
      break;
    case  2:  // SOUTH        y  
      GlobalBC::bc_fixP(0.0,1.0,0.0,&closures,prob_parm.p_0, s_int, s_ext); 
      break;
    case -2:  // NORTH        y 
      GlobalBC::bc_fixP(0.0,-1.0,0.0,&closures,prob_parm.p_0, s_int, s_ext); 
      break;
    default:

      break; 
  }
}
///////////////////////////////AMR//////////////////////////////////////////////
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
user_tagging(int i, int j, int k, int nt_level, auto &tagfab,
             const auto &sdatafab, const auto& ebflag, const auto &geomdata,
             const ProbParm &prob_parm, int level) {

  const Real *prob_lo = geomdata.ProbLo();
  const Real *dx = geomdata.CellSize();
  Real x = prob_lo[0] + (i + Real(0.5)) * dx[0];
  Real y = prob_lo[1] + (j + Real(0.5)) * dx[1];
  Real z = prob_lo[2] + (k + Real(0.5)) * dx[2];
  Real r = sqrt(x*x + y*y);

  bool refine = false;
        
  switch (level)

    {
      case 0:
        refine = (x < 1.4e-3) && (y < 1.1e-3);
        //refine = true;
        break;
      case 1:
        refine = (x > 0) && (y < 1.038e-3) && (x < 0.0012972);

        break;
      case 2:
        refine= (x > 4.350e-4) && (x < 6.690e-4) && (y < 0);


        break;
      case 3:
        // refine= (z > 0.035) && (z < 0.07);    
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
  bool static constexpr ATF = false; // use adaptive thickening factor
  static constexpr Real thickfactor = 0.0; // thickening factor

  bool static constexpr do_reactions = false;

  // viscous options
  static constexpr int order = 2;                  // order numerical scheme   
  static constexpr bool use_LES= true;
  bool static constexpr mask_cells_boundary = true; // avoid compute reactions in cells partially covered

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
  void inline src(const Geometry& geomdata, const amrex::MFIter &mfi,
                  const amrex::Array4<const amrex::Real> &prims,
                  const amrex::Array4<amrex::Real> &rhs, const cls_t *cls_d,
                  amrex::Real dt, amrex::Real real_time, const auto& ebflag){

    const Box& bxg = mfi.tilebox();
    // const Box& bxg = mfi.growntilebox(cls_t::NGHOST);
    const Real *prob_lo = geomdata.ProbLo();
    const Real *dx = geomdata.CellSize();

    ProbParm const prob_parm;
    const auto& cls = *cls_d;

    // buffer coordinates.
    // west wall (inflow)
    const Real x0_lo = 7.36e-5; 
    const Real xw_lo = -0.0001091;
    // east wall (outflow)
    const Real x0_hi = 0.0011145;
    const Real xw_hi = 0.0012972;
    // top wall (normal)
    const Real y0_hi = 0.998e-3;
    const Real yw_hi = 0.0011807;

    amrex::ParallelFor(bxg,
      [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
      
      const Real x = prob_lo[0] + (i + Real(0.5)) * dx[0];
      const Real y = prob_lo[1] + (j + Real(0.5)) * dx[1];
      const Real z = prob_lo[2] + (k + Real(0.5)) * dx[2];

      bool buffer = (x <= x0_lo) || (x >= x0_hi) || (y >= y0_hi);

      if (buffer) {
        // 1. Calculate normalized penetration depth (0 to 1) for each boundary
        // Note: Assumes xw_hi > x0_hi and yw_hi > y0_hi
        Real d_xlo = (x <= x0_lo) ? (x0_lo - x) / (x0_lo - xw_lo) : 0.0;
        Real d_xhi = (x >= x0_hi) ? (x - x0_hi) / (xw_hi - x0_hi) : 0.0;
        Real d_yhi = (y >= y0_hi) ? (y - y0_hi) / (yw_hi - y0_hi) : 0.0;

        // 2. Combine using max to naturally handle corners
        Real a_linear = amrex::max(d_xlo, amrex::max(d_xhi, d_yhi));
        a_linear = amrex::min(1.0, a_linear);
        Real tau_relax = a_linear * 200000; //corresponds to tau relax = 5.0e-5 and coef = dt/taurelax
        Real coef = dt * tau_relax;

        // pressure relax  if P > P0  P drops and to keep T constant rho drops
        Real pres = prims(i,j,k,cls.QPRES);

        const Real dP = (prob_parm.p_0- pres)*coef;
      
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
        rhs(i,j,k,cls.URHO) += drhodt;     
        rhs(i,j,k,cls.UMX) += prims(i,j,k,cls.QU)*drhodt;
        rhs(i,j,k,cls.UMY) += prims(i,j,k,cls.QV)*drhodt; 
        rhs(i,j,k,cls.UMZ) += prims(i,j,k,cls.QW)*drhodt;
        rhs(i,j,k,cls.UET) += Et*drhodt; 
      }
    });

  };
};

} // namespace PROB

#endif
