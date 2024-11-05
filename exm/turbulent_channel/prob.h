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

// Periodic Turbulent channel

using namespace amrex;


namespace PROB {


static constexpr Real Reynolds = 3000.0;  // bulk Reynolds number
static constexpr Real Mach     = 1.5;     // bulk Mach number
static constexpr Real Rgas     = gas_constant*1000.0/28.96;  // bulk Mach number
static constexpr Real Ggas     = 1.4;      // gamma
static constexpr Real viscos   = 1.0/Reynolds;    // viscos temp
static constexpr Real Retau    = 100;      // Retau



// problem parameters 
struct ProbParm {
  Real Tw   =   500;	
  Real h    =   0.006845;     // channel height 2 h
  Real ubulk    = Mach*sqrt(Ggas*Rgas*Tw); // bulk velocity
  //Real mass     = Reynolds*viscos/h;       // mass flow rate
  Real rhob     = 1.0;      
  Real Q        = ubulk*8*h*h;
  Real mass     = rhob*ubulk;
  Real utau     = ubulk*Retau/Reynolds;
  Real tau      = rhob*utau*utau;          // tau wall
  Real p        = rhob*Rgas*Tw;            // pressure        
  Real L        = 2*h;
  // Real dpdx     = 0.5*tau/(rhob*h);       // forcing term (pressure gradient)
  Real dpdx = 12.0*Q*Q/(Reynolds*L*L*L);    // forcing term (pressure gradient)
};

// numerical method parameters
struct methodparm_t {

  public:

  static constexpr int  order = 2;                   // order numerical scheme viscous
  static constexpr Real conductivity = 0.0262;       // conductivity (for constant value)
  static constexpr Real viscosity    = viscos;        // viscosity    (for constant value)

  static constexpr bool solve_viscterms_only = false;
  
};

struct skewparm_t {

  public:

  static constexpr int  order = 4;                 
  static constexpr bool dissipation = true;  
  static constexpr bool ibm = false;
  static constexpr Real C2skew=0.1,C4skew=0.0016; 
  
};


inline Vector<std::string> cons_vars_names={"Xmom","Ymom","Zmom","Energy","Density"};
inline Vector<int> cons_vars_type={1,2,3,0,0};

//typedef closures_dt<indicies_t, visc_const_t<methodparm_t>, cond_const_t<methodparm_t>,
//                    calorifically_perfect_gas_t<indicies_t>> ProbClosures;

typedef closures_dt<indicies_t, transport_const_t<methodparm_t>,
                    calorifically_perfect_gas_t<indicies_t>> ProbClosures;


template <typename cls_t > class user_source_t;

typedef rhs_dt<skew_t<skewparm_t,ProbClosures>, viscous_t<methodparm_t, ProbClosures>, user_source_t<ProbClosures> > ProbRHS;



void inline inputs() {
  ParmParse pp;

  pp.add("cns.order_rk", 3);   // -2, 1, 2 or 3"
  pp.add("cns.stages_rk", 3);  // 1, 2 or 3

  amrex::Print() << " ****** Starting ... *******" <<  std::endl;
  amrex::Print() << " 3D Turbulent Channel Flow  (Oct 2024)" <<  std::endl;

 // std::mt19937 generator(rd());

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

  const Real yoh = y/prob_parm.L;
  const Real Lx  = 12.0*prob_parm.h;  
  const Real Lz  = 4.0*prob_parm.h;

  
  // local vars
  Real rho,eint, T, P,u[3],ufluc[3]={0.0, 0.0, 0.0};  
  
  std::random_device rd;
  std::mt19937 generator(rd());
  std::uniform_real_distribution<double> distribution(-1.0,1.0);
  std::uniform_real_distribution<double> distribution2(-0.01,0.01);
  double rn = distribution(generator);

  // fluctuations
  
  // initial conditioms  
  P    = prob_parm.p; rho = prob_parm.rhob;
  Real Ub = prob_parm.ubulk; 
  u[1] = Real(0.0); u[2] = Real(0.0);
  u[0] = 12.0*Ub*(1.0 - yoh)*yoh;
  
  // fluctuations
  const int NFREQ=6;
  Real freqz[NFREQ],Amp[NFREQ],freqx[NFREQ];
  for (int l=0;l<NFREQ;l++) 
  {
    double r1 = distribution2(generator);
    double r2 = distribution2(generator);
    //r1= 0.0; r2=0.0;
    Amp[l]   = 0.05/(l+1);
    freqx[l] = l*2.0*M_PI*x/Lx + r1*M_PI;
    freqz[l] = l*2.0*M_PI*z/Lz + r2*M_PI;
  }

  // fluctuations
  ufluc[0]+= 0.4*Ub*rn;
  for (int l=0;l<NFREQ;l++)
  {
    ufluc[1]+= Amp[l]*Ub*sin(freqx[l])*cos(freqz[l] );
    ufluc[2]+= Amp[l]*Ub*cos(freqx[l])*sin(freqz[l] );
  }


  //printf(" z=%f vf=%f wf=%f \n ",z, ufluc[1],ufluc[2]);
  
  Real kin =0.0;
  for (int dim=0;dim<3;dim++)
  {u[dim]+=ufluc[dim]; kin+= Real(0.5)*u[dim]*u[dim];}


  // T and E from P and T
  T =   P/(cls.Rspec*rho);
  eint =  cls.cv*T;

  state(i, j, k, cls.URHO) = rho;
  state(i, j, k, cls.UMX)  = rho * u[0];
  state(i, j, k, cls.UMY)  = rho * u[1];
  state(i, j, k, cls.UMZ)  = rho * u[2];  
  state(i, j, k, cls.UET)  = rho*eint + kin;    
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

  // rho, e ---> P
  Real rho  = s_int[URHO];
  Real rhoinv = Real(1.0) / rho;
  Real ux = s_int[UMX] * rhoinv;
  Real uy = s_int[UMY] * rhoinv;
  Real uz = s_int[UMZ] * rhoinv;
  Real rhoke = Real(0.5) * rho * (ux * ux + uy * uy + uz * uz);
  Real rhoei = (s_int[UET] - rhoke);
  Real p = (Ggas - Real(1.0)) * rhoei;
  // comoute rhow
  Real rhow = p/(Rgas*prob_parm.Tw);
  Real rho1 = 2.0*rhow-rho;

  

  switch(face)
  {
    case 2: // SOUTH adiabatic wall
      s_ext[URHO] = rho1;
      s_ext[UMX]  = -s_int[UMX];
      s_ext[UMY]  = -s_int[UMY];
      s_ext[UMZ]  = -s_int[UMZ];
      s_ext[UET]  = s_int[UET];
      break;
    case -2: // NORTH adiabatic wall
      s_ext[URHO] = rho1;
      s_ext[UMX]  = -s_int[UMX];
      s_ext[UMY]  = -s_int[UMY];
      s_ext[UMZ]  = -s_int[UMZ];
      s_ext[UET]  = s_int[UET];
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

      //  Pressure Source           
      //Real rho = prims(i, j, k, cls.QRHO);
      rhs(i,j,k,cls.UMX) += prob_parm.dpdx;  

      });

  };
};

///////////////////////////////AMR//////////////////////////////////////////////
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
user_tagging(int i, int j, int k, int nt_level, auto &tagfab,
             const auto &sdatafab, const auto &geomdata,
             const ProbParm &prob_parm, int level) {


  Real Threshold[3] = {100.0,100.0,100.0};
  
  Real rhop = sdatafab(i,j,k,ProbClosures::URHO);
  Real drhox = Math::abs(sdatafab(i+1,j,k,ProbClosures::URHO) -
              sdatafab(i-1,j,k,ProbClosures::URHO));
  Real drhoy = Math::abs(sdatafab(i,j+1,k,ProbClosures::URHO) -
              sdatafab(i,j-1,k,ProbClosures::URHO));
  Real gradrho= sqrt(drhox*drhox+drhoy*drhoy)/rhop;

  if (nt_level> 0)
  {
    //tag cells based on density  values       
    tagfab(i,j,k) = (gradrho > Threshold[level]);        
      
  }
}
////////////////////////////////////////////////////////////////////////////////




///////////////////////////////////////////////////////////////////////////////

} // namespace PROB
#endif
