// euler equations except at walls

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

#include <bc_types.h>


// NTU Combustor-type for demonstration purposes
// created by S Dupre and S Navarro-Martinez (2025)


using namespace amrex;

namespace PROB {

static constexpr Real Mw   = 28.96e-3;    // Molecular weight
static constexpr Real gam  = 1.4;         // Adiabatic coefficient

// problem parameters 

struct ProbParm {  
  // inflow state
  static constexpr Real p_inflow    = pres_atm2si; //[Pa] inflow pressure (1 atm)  
  static constexpr Real T_inflow    = 285.5;  //[K]  
  static constexpr Real rho_inflow  = p_inflow*Mw/(gas_constant*T_inflow);
  static constexpr Real vel_in[3]= {0.0,0.0,7.9}; // array of inflow velocity [m/s]
  static constexpr Real eint_inflow = p_inflow/ (gam - Real(1.0));  

  static constexpr Real rhou_inflow = rho_inflow*vel_in[0];
  static constexpr Real rhov_inflow = rho_inflow*vel_in[1];
  static constexpr Real rhow_inflow = rho_inflow*vel_in[2];
  static constexpr Real kin = Real(0.5) * rho_inflow *
      (vel_in[0]*vel_in[0]+vel_in[1]*vel_in[1]+vel_in[2]*vel_in[2]);
  static constexpr Real rhoe_inflow = eint_inflow + kin;  
  
  // inside combustor state/exit
  static constexpr Real p_0     = pres_atm2si; //[Pa] inflow pressure (1 atm) 
  static constexpr Real T_0     = 285.5;  //[K]  
  static constexpr Real rho_0   = p_0*Mw/(gas_constant*T_0);
  static constexpr Real vel_0[3]= {0.0,0.0,0.0}; // array of inside velocity [m/s]
  static constexpr Real eint_0  = p_0/ (gam - Real(1.0));  

  // volumetric flow rate
  static constexpr Real Q = 9.35;  // rho U   kg /m2 s

  // species array
  static constexpr Real Y_inflow[NUM_SPECIES] = {1.0};


  // geometry auxiliary
  Real zin = 0.01;
  Real zexit = 0.135;
  
};


// numerical method parameters
struct const_viscparm_t {

  public:

  //static constexpr bool dissipation = true;         // no dissipation
  //static constexpr int  order = 4;                  // order numerical scheme   
  //static constexpr Real C2skew=0.5,C4skew=0.016;   // Skew symmetric default
  static constexpr Real viscosity = 1.846e-5;
  static constexpr Real conductivity = 0.02624;
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


typedef closures_dt<indicies_t, visc_suth_t, cond_suth_t, calorifically_perfect_gas_t<indicies_t>> ProbClosures;

template <typename cls_t > class user_source_t;

//typedef closures_dt<indicies_t, transport_const_t<const_viscparm_t>, calorifically_perfect_gas_t<indicies_t>> ProbClosures;

// define nuemrical scheme comment/uncomment to set up 
typedef rhs_dt<weno_t<ReconScheme::Teno5, ProbClosures>, viscous_t<viscous_param_t, ProbClosures>, user_source_t<ProbClosures> > ProbRHS;
//typedef rhs_dt<riemann_t<false, ProbClosures>, viscous_t<viscous_param_t, ProbClosures>, user_source_t<ProbClosures> > ProbRHS;


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

  Real rad = sqrt(x*x + y*y);

  // local vars
  Real rhot,eint,u[3];

  
  rhot =  prob_parm.rho_inflow;    
  for(int idim=0;idim < AMREX_SPACEDIM;idim++) {u[idim]=prob_parm.vel_in[idim];}
  eint =  prob_parm.eint_inflow;

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
      GlobalBC::bc_fixP(0.0,1.0,0.0,&closures,prob_parm.p_0, s_int, s_ext); 
      break;
    case  1:  // WEST
      GlobalBC::bc_fixP(1.0,0.0,0.0,&closures,prob_parm.p_0, s_int, s_ext); 
      break;
    case -1:  // EAST
      GlobalBC::bc_fixP(-1.0,0.0,0.0,&closures,prob_parm.p_0, s_int, s_ext);  
      break;
    case -2:  // NORTH
      GlobalBC::bc_fixP(-1.0,0.0,0.0,&closures,prob_parm.p_0, s_int, s_ext); 
      break;
    case -3:   //RIGHT 
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
  const Real r2 = sqrt(x*x + y*y + (z-prob_parm.zexit)*(z-prob_parm.zexit));
  


  bool refine = false;
  //std::cout << "refined!" << std::endl; 
  //refine = (z < prob_parm.zexit); 


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
      std::cout << "refined!" << std::endl; 
      refine = (z < prob_parm.zexit);
      //refine = (((z < prob_parm.zexit) || (r2 < 0.025) ) && (r < 0.025) )  ;    // refine combustor    
      break;
    case 1:
      std::cout << " second level reached " << std::endl;
      //refine= (z > 0.035) && (z < 0.07);
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
  void inline src(const Geometry& geomdata, const amrex::MFIter &mfi,
                  const amrex::Array4<const amrex::Real> &prims,
                  const amrex::Array4<amrex::Real> &rhs, const cls_t *cls_d,
                  amrex::Real dt){

    //const Box bx = mfi.tilebox();
    const Box& bxg = mfi.growntilebox(cls_t::NGHOST);
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
      const Real tau_relax = 50.0*dt; const Real coef =dt/tau_relax;

      const Real dP = (prob_parm.p_0- prims(i,j,k,cls.QPRES))*coef;
      
      const Real T = prims(i,j,k,cls.QT); // dT =0
      const Real rho  = prims(i, j, k, cls.QRHO);
      Real drho = 0.0; const Real Y[NUM_SPECIES] = {1.0};
      cls.PYT2R(dP,Y,T,drho);  Real drhodt = drho/dt;
      Real kin = 0.5*(prims(i,j,k,cls.QU)*prims(i,j,k,cls.QU) + prims(i,j,k,cls.QV)*prims(i,j,k,cls.QV)
                    + prims(i,j,k,cls.QW)*prims(i,j,k,cls.QW));
      Real Et  = prims(i,j,k,cls.QEINT) + kin;

      bool buffer = (z > prob_parm.zexit) && (r > 0.035);

      if (buffer){        
        rhs(i,j,k,cls.URHO)+= drhodt;
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
