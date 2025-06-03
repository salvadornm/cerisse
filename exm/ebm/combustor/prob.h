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


  // geometry auxiliary
  Real zin = 0.01;
};


// numerical method parameters
struct const_viscparm_t {

  public:

  //static constexpr bool dissipation = true;         // no dissipation
  //static constexpr int  order = 5;                  // order numerical scheme   
  //static constexpr Real C2skew=0.5,C4skew=0.0016;   // Skew symmetric default
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
      // parameters from interior flow
	    {
      
      Real q_int[ProbClosures::NPRIM]={0.0};      
      closures.cons2prims_point(s_int,q_int); // convert to primitive

      // pghost =pinterior      
      const Real P = q_int[ProbClosures::QPRES];            
      const Real rho  = s_int[URHO] * q_int[ProbClosures::QT]/ prob_parm.T_inflow;  
      s_ext[URHO] = rho;

      // ensure mass flow rate is constant - no acoustic reflection at inlet      
      s_ext[UMX] = 0;
      s_ext[UMY] = 0;
      s_ext[UMZ] = prob_parm.Q;

      // total energy per unit vol
      Real e_ext = 0.0;
      const Real Y[NUM_SPECIES] = {1.0};
      closures.RYP2E(rho, Y, P, e_ext);      
      s_ext[UET] = rho* e_ext + 0.5 * s_ext[UMZ] * s_ext[UMZ] /rho;
      
      /**  
      s_ext[URHO] = prob_parm.rho_inflow;
      s_ext[UMX]  = prob_parm.rhou_inflow;
      s_ext[UMY]  = prob_parm.rhov_inflow;
      s_ext[UMZ]  = prob_parm.rhow_inflow;
      s_ext[UET]  = prob_parm.rhoe_inflow;
      */
      break;
      }
    case  2:  // SOUTH
      break;
    case  1:  // WEST
      break;
    case -1:  // EAST
      break;
    case -2:  // NORTH
      break;
    
    case -3:   //RIGHT 
      {
      Real q_int[ProbClosures::NPRIM]={0.0};      
      closures.cons2prims_point(s_int,q_int); // convert to primitive
      // Pout -> P0  
      const Real P = prob_parm.p_0; 
    
      // Tout = Tint  (dT/dz =0)
      const Real T = q_int[ProbClosures::QT];
      
      Real Y[NUM_SPECIES] = {1.0};
      Real rho = 0.0; Real e_ext=0.0;
      closures.PYT2R(P,Y,T,rho);
      closures.PYT2E(P,Y,T,e_ext);

      s_ext[URHO] = rho;
      s_ext[UMX]  = rho*q_int[ProbClosures::QU];
      s_ext[UMY]  = rho*q_int[ProbClosures::QV];
      s_ext[UMZ]  = rho*max(q_int[ProbClosures::QW],0.0);

      const Real kin = 0.5*(s_ext[UMX] * s_ext[UMX] + s_ext[UMY] * s_ext[UMY] + s_ext[UMZ] * s_ext[UMZ]);

      s_ext[UET] = rho* e_ext + kin/rho;

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

  // refine exit of injector
  bool refine= (z > 0.035) && (z < 0.07);

  // refine close to exit  (avoid corner problem)
  refine= (z > 0.13) || refine;


  //const int URHO= ProbClosures::URHO;

  // // compoute | d rho | normalised with rho
  // Real o_over_rhot = Real(1.0)/sdatafab(i,j,k,URHO);

  // Real drhox = Math::abs(sdatafab(i+1,j,k,URHO) - sdatafab(i-1,j,k,URHO));
  // Real drhoy = Math::abs(sdatafab(i,j+1,k,URHO) - sdatafab(i,j-1,k,URHO));
  // Real drhoz = Math::abs(sdatafab(i,j,k+1,URHO) - sdatafab(i,j,k-1,URHO));

  // Real gradrho= Real(0.5)*sqrt(drhox*drhox+drhoy*drhoy)*o_over_rhot;        


  // switch (level)
  // {
  //   case 0:
  //     tagfab(i,j,k) = (gradrho > 0.1);        
  //     break;
  //   case 1:
  //     tagfab(i,j,k) = (gradrho > 0.2);        
  //     break;
  //   default:
  //     tagfab(i,j,k) = (gradrho > 0.3);        
  //   break;
  //  }
    

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
      
      const Real z    = prob_lo[2] + (k + Real(0.5)) * dx[2];
      
      const Real Lout =  0.5;
      const Real fz =  (prims(i,j,k,cls.QPRES) - prob_parm.p_0)/Lout;  // dp/dx
      // if P > P0 small aceleration

      // pressure relax  if P > P0  P drops and to keep T constant rho drops
      const Real tau_relax = 25.0*dt;
      const Real dP = (prob_parm.p_0- prims(i,j,k,cls.QPRES))*dt/tau_relax;

      //const Real Pnew = prims(i,j,k,cls.QPRES) + dP;
      const Real Tnew = prims(i,j,k,cls.QT); // dT =0
      const Real rho  = prims(i, j, k, cls.QRHO);
      Real drho = 0.0; Real Y[NUM_SPECIES] = {1.0};
      cls.PYT2R(dP,Y,Tnew,drho); 
      Real drhodt = drho/dt;
      const Real kin = prims(i,j,k,cls.QU)*prims(i,j,k,cls.QU) + 
      prims(i,j,k,cls.QV)*prims(i,j,k,cls.QV) + prims(i,j,k,cls.QW)*prims(i,j,k,cls.QW);
      const Real Et  = prims(i,j,k,cls.QEINT) + 0.5*kin;
      Real edesired = 0.0;
      cls.PYT2E(prob_parm.p_0,Y,300.0,edesired);
      if (z > 0.125){        
        rhs(i,j,k,cls.URHO)+= drhodt;
        //rhs(i,j,k,cls.UMX) += prims(i,j,k,cls.QU)*drhodt;
        //rhs(i,j,k,cls.UMY) += prims(i,j,k,cls.QV)*drhodt;
        rhs(i,j,k,cls.UMZ) += rho*fz;// + prims(i,j,k,cls.QW)*drhodt;
        //rhs(i,j,k,cls.UET) += rho*fz*prims(i,j,k,cls.QW);//
        rhs(i,j,k,cls.UET) += Et*drhodt;

        rhs(i,j,k,cls.UET) += (edesired - prims(i,j,k,cls.QEINT))*dt/tau_relax;

      }


      });

  };
};


} // namespace PROB

#endif
