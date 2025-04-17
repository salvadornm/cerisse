#ifndef CNS_PROB_H_
#define CNS_PROB_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_ParmParse.H>
#include <Closures.h>
#include <RHS.h>

#include <Constants.h>
#include <NumParam.h>

#include <random>

#if CNS_USE_EB    
#include <ebm.h>
#include <walltypes.h>
#endif

// Turbulent Pipe Flow
// created by S Navarro-Martinez (2025)

using namespace amrex;

namespace PROB {

static constexpr Real Mw       = 28.96e-3;         // Molecular weight
static constexpr Real gam      = 1.4;              // Adiabatic coefficient
static constexpr Real Rgas     = gas_constant/Mw;  // bulk Mach number
static constexpr Real Reynolds = 24600;            // bulk Reynolds number
static constexpr Real Pr       = 0.7;              // Prandtl
static constexpr Real Cv       = Rgas/(gam - 1.0);
static constexpr Real Cp       = gam*Cv;
static constexpr Real viscos   = 1.0/Reynolds;    // constant viscosity
static constexpr Real lambda   = viscos*Cp/Pr;    // constant lambda
static constexpr Real Retau    = 1920;            // (based on friction velocity)


// problem parameters 

struct ProbParm {  
  
  // inside pipe
  Real p_0     = pres_atm2si; //[Pa] inflow pressure (1 atm) 
  Real T_0     = 300.0;  //[K]  
  Real rho_0   = p_0/(Rgas*T_0);
  Real vel_0[3]= {0.0,0.0,1.0}; // array of inside velocity [m/s]
  Real eint_0  = p_0/ (gam - Real(1.0));  

  // geometry auxiliary
  Real zin = 0.01;
  Real Dpipe = 1.0;
  Real Rpipe = 0.5*Dpipe;

  // forcing
  Real tau_w = 0.00284;
  Real utau = sqrt(tau_w/rho_0);
  Real fx = tau_w/(Rpipe*rho_0);       // forcing term per unit mass 

};

// numerical method parameters
struct methodparm_t {

  public:

  static constexpr int  order = 2;                  // order numerical scheme   
  static constexpr Real conductivity = lambda;       // conductivity (for constant value)
  static constexpr Real viscosity    = viscos;       // viscosity    (for constant value)

};

// parmeters for skew-symmetric method
struct skewparm_t {

  public:

  static constexpr bool dissipation = true;         // no dissipation
  static constexpr int  order = 4;                  // order numerical scheme   
  static constexpr Real C2skew=0.5,C4skew=0.0016;   // Skew symmetric default

};

inline Vector<std::string> cons_vars_names={"Xmom","Ymom","Zmom","Energy","Density"};
inline Vector<int> cons_vars_type={1,2,3,0,0};

// closures
typedef closures_dt<indicies_t, transport_const_t<methodparm_t>,calorifically_perfect_gas_t<indicies_t> > ProbClosures;

template <typename cls_t > class user_source_t;

// building rhs
//typedef rhs_dt<skew_t<skewparm_t,ProbClosures>, no_diffusive_t, user_source_t<ProbClosures> > ProbRHS;

typedef rhs_dt<skew_t<skewparm_t,ProbClosures>, viscous_t<methodparm_t, ProbClosures>, user_source_t<ProbClosures> > ProbRHS;
//typedef rhs_dt<weno_t<ReconScheme::Teno5, ProbClosures>, no_diffusive_t, no_source_t > ProbRHS;


// define type of wall and EBM class
#if CNS_USE_EB    
typedef adiabatic_wall_t<ProbClosures> TypeWall;
typedef ebm_t<TypeWall,ProbClosures> ProbEB;
#endif

void inline inputs() {

  amrex::Print() << " ****** Starting ... ******* " <<  std::endl;
  amrex::Print() << " 3D Turbulent Pipe Flow  (Apr 2025)" <<  std::endl;
  
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
  Real r =  rad/prob_parm.Rpipe;

  // local vars
  Real rhot,eint,u[3],ufluc[3]={0.0, 0.0, 0.0};

  rhot =  prob_parm.rho_0;    
  for(int idim=0;idim < AMREX_SPACEDIM;idim++) {u[idim]=prob_parm.vel_0[idim];}
  eint =  prob_parm.eint_0;

  if (r < 1) 
  {
    u[2] = prob_parm.vel_0[2]*(1.0-r*r);
  }

  // fluctuations
  std::random_device rd;       
  std::mt19937 generator(rd());
  std::uniform_real_distribution<double> distribution(-1.0,1.0);
  // 
  const Real Lz = 4.8;
  const int NFREQ=6;
  Real freqz[NFREQ],Amp[NFREQ];
  for (int l=0;l<NFREQ;l++)
  {
    double r1 = 0.01*distribution(generator); //random shift
    Amp[l]   = 0.1/(l+1);
    freqz[l] = l*2.0*M_PI*x/Lz + r1*M_PI;
    ufluc[0] +=Amp[l]*sin(freqz[l]);
    ufluc[1] +=Amp[l]*cos(freqz[l]); 
  }
  double rn = distribution(generator);
  ufluc[2]+= 0.02*rn;  

  //
  // add fluctuations to mean flow
  for (int dim=0;dim<3;dim++)
  {u[dim]+=ufluc[dim];}



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
      break;
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
///////////////////////////////SOURCE TERM /////////////////////////////////////
template <typename cls_t>
class user_source_t {
  public:
  void inline src(const amrex::MFIter &mfi,
                  const amrex::Array4<const amrex::Real> &prims,
                  const amrex::Array4<amrex::Real> &rhs, const cls_t *cls_d,
                  amrex::Real dt){

    const Box bx = mfi.tilebox();

    ProbParm const prob_parm;

    amrex::ParallelFor(bx,
      [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
      const auto& cls = *cls_d;

      Real rho = prims(i, j, k, cls.QRHO);
      rhs(i,j,k,cls.UMX) += rho*prob_parm.fx;
      rhs(i,j,k,cls.UET) += rho*prob_parm.fx*prims(i,j,k,cls.QU);

      });

  };
};
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
  Real rad = sqrt(x*x + y*y);
  Real r =  rad/prob_parm.Rpipe;

  // const int URHO= ProbClosures::URHO;

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
    
  // refine close to pipe

  tagfab(i,j,k) = (r > 0.9);        

}
///////////////////////////////////////////////////////////////////////////////

} // namespace PROB

#endif
