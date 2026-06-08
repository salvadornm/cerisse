#ifndef CNS_PROB_H_
#define CNS_PROB_H_

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_ParmParse.H>
#include <PelePhysics.H>
#include <ReactorBase.H> 

#include <Closures.h>
#include <RHS.h>
#include <ibm_solver.h>
#include <ratio>
#include <Constants.h>
#include <ibm_walltypes.h>

#include <ibm_solver.h>
#include <Constants.h>
#include <NozzleFunctions.h>
#include <ibm_wallmodel.h>
#include <numbers>

using namespace amrex;
using namespace universal_constants;

namespace PROB {

// constants  
static constexpr Real Mach     = 2.0;              // bulk Mach number
static constexpr Real Mw       = 44.e-3;           // Molecular weight [CO2]
static constexpr Real gam      = 1.3;              // Adiabatic coefficient 
static constexpr Real gam_srp  = 1.4;              // Adiabatic coefficient 
static constexpr Real Rgas     = gas_constant/Mw;  // gas constant
static constexpr Real Reynolds = 10000;            // bulk Reynolds number
static constexpr Real Pr       = 0.7;              // Prandtl
static constexpr Real Cv       = Rgas/(gam - 1.0);
static constexpr Real Cp       = gam*Cv;
static constexpr Real viscos   = 1.0/Reynolds;    // constant viscosity
static constexpr Real lambda   = viscos*Cp/Pr;    // constant lambda

static constexpr int ibm_eorder=1;
  
//////////////////////////// Physical modelling ////////////////////////////////
struct ProbParm
{ 

  //------ values obtained from ./properties.py  (SI units)---------
  // based on Mach number 2.0 and gamma 1.334533782996543
  // and free-stream conditions P_oo = 284.0 [Pa] and T_oo = 227.0 [K]
  // Staganation conditions
  static constexpr Real P0  =  2191.8913255191455 ;
  static constexpr Real T0  =  378.8783374804305 ;
  // Free-stream conditions
  static constexpr Real p_oo  =  284.0 ;
  static constexpr Real T_oo  =  227.0 ;
  static constexpr Real rho_oo  =  0.00659535789648574 ;
  static constexpr Real c_oo  =  239.7202778419921 ;
  static constexpr Real u_oo  =  479.4405556839842 ;
  static constexpr Real eint_oo  =  -8683464.706718355 ;
  static constexpr Real kin_oo  =  758.0151887420365 ;
  // SRP Staganation conditions
  static constexpr Real P0srp  =  2191.8913255191455 ;
  static constexpr Real T0srp  =  1894.3916874021525 ;
  // SRP throat conditions (choked flow)
  static constexpr Real Pt  =  1195.1513931292714 ;
  static constexpr Real Tt  =  1645.445074309705 ;
  static constexpr Real u_srp  =  855.7917341280771 ;
  //-----------------------------------------------------------------
  
  static constexpr Real  Yco2_oo = 0.96; // 96% CO2
  static constexpr Real  Yar_oo  = 0.04; //  4% Argon
  
  // centre of probe (approx)
  static constexpr Real x0 = 1.0, y0 = 2.0, z0 = 2.0;
  static constexpr Real xsrp= 0.6225,ysrp = 0.75, zsrp=0.75;
  static constexpr Real Rsrp= 0.008; // nozzle radius
  static constexpr Real Asrp= std::numbers::pi*Rsrp*Rsrp;  
};

//  parameters for viscous solver and conductivity/viscosity
struct methodparm_t {

  public:

  static constexpr int  order = 2;                  // order numerical scheme   
  static constexpr Real conductivity = lambda;       // conductivity (for constant value)
  static constexpr Real viscosity    = viscos;       // viscosity    (for constant value)
};


// parameters for skew-symmetric method
struct skewparm_t {

  public:

  static constexpr bool dissipation = true;         // no dissipation
  static constexpr int  order = 4;                  // order numerical scheme   
  static constexpr Real C2skew=1.5,C4skew=0.016;   // Skew symmetric values 
};

struct ibmparm_t {

  public:

  static constexpr int  interp_order = 1;
  static constexpr int  extrap_order = ibm_eorder;
  static constexpr Real alpha= 0.6; 
  // surface parameters
  static constexpr int  interp_order_surf = 1;
  static constexpr int  extrap_order_surf = 1;
  static constexpr Real alpha_surf= 0.6;     
  // 
  static constexpr int  ghost_layers = 1;
  static constexpr bool interior_is_solid = true;
};

// CLOSURES
using ProbClosures = closures_dt< indicies_t, transport_Pele_t, multispecies_pele_gas_t<indicies_t> >;


// NUMERICAL SCHEME + EQNS TO SOLVE   (Euler/NS/Source)                 

typedef rhs_dt<riemann_t<false, ProbClosures>, no_diffusive_t, no_source_t > ProbRHS;

//using ProbRHS = rhs_dt< riemann_t<false, ProbClosures>, viscous_t<methodparm_t, ProbClosures>, reactor_t<ProbClosures> >;

// declaration of user-specific ibm (see bottom of the file for definition)
template < typename param, typename cls_t > class ibm_user_t; 

// IBM templates
typedef ibm_user_t<ProbParm,ProbClosures> TypeWall;
typedef ibm_solver_t<TypeWall,ibmparm_t,ProbClosures> ProbIB;


void inline inputs() {
  
  amrex::Print() << " ****** Starting ... ******* " <<  std::endl;
  amrex::Print() << " Supersonic Flow over Phoebus geom " <<  std::endl;
  amrex::Print() << " free stream mimic mars atmosphere " <<  std::endl;
  amrex::Print() << " SRP jet of Nitrogen " <<  std::endl;  
  amrex::Print() << " **************************** " <<  std::endl;
}

//////////////////////////// Initial conditions ////////////////////////////////
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void prob_initdata (int i, int j, int k, amrex::Array4<amrex::Real> const& state,
      amrex::GeometryData const& geomdata, ProbClosures const& cls, ProbParm const& pparm) {
  
  //const Real* prob_lo = geomdata.ProbLo();  
  //const Real* dx      = geomdata.CellSize();

  // Real x = prob_lo[0] + (i+0.5_rt)*dx[0];
  // Real y = prob_lo[1] + (j+0.5_rt)*dx[1];
  // Real z = prob_lo[2] + (k+0.5_rt)*dx[2];

  // local vars
  Real rhot,eint,u[3]={0.0};
  Real Yt[NUM_SPECIES] ={0.0};

  // initial state
  rhot =  pparm.rho_oo;
  u[0] =  pparm.u_oo;
  eint =  pparm.rho_oo*pparm.eint_oo;
  
  Yt[CO2_ID]   = pparm.Yco2_oo;
  Yt[AR_ID]    = pparm.Yar_oo;

  state(i, j, k, cls.UMX)  = rhot * u[0];
  state(i, j, k, cls.UMY)  = rhot * u[1];
  state(i, j, k, cls.UMZ)  = rhot * u[2];
  state(i, j, k, cls.UET)  = eint + Real(0.5) * rhot * u[0] * u[0] ; 
  for (int n = 0; n < NUM_SPECIES; ++n) {
    state(i, j, k, cls.UFS + n) = rhot * Yt[n];
  }
  

}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE 
void user_tagging(int i, int j, int k, int nt, auto& tagfab, const auto &sdatafab, 
                  const Array4<const unsigned char>& ibfab, const auto& geomdata, 
                  const ProbParm& pparm , int level) {

  const Real* dx  = geomdata.CellSize();
  const Real x = (i+0.5_rt)*dx[0];
  const Real y = (j+0.5_rt)*dx[1];
  const Real z = (k+0.5_rt)*dx[2];
  Real xrel[3];
  // coordinate relative to object
  xrel[0]= x-pparm.x0; xrel[1]= y-pparm.y0; xrel[2]= z-pparm.z0;
  Real radius =xrel[0]*xrel[0] + xrel[1]*xrel[1] + xrel[2]*xrel[2];
  const Real Rmax = 0.5_rt*0.5_rt;const Real Rmin = 0.15_rt*0.15_rt;
  // initialize thresholds at all levels
  Real rhofluc_threshold[6] = {0.3_rt,0.6_rt,0.9_rt,1000_rt,1000_rt,1000_rt};

  
  int URHO = ProbClosures::UFS + 1; // rho ~ rho*YCO2 
  // refine close to grads of rho*YCO2 
  Real drhox = std::abs(sdatafab(i+1,j,k,URHO) - sdatafab(i-1,j,k,URHO))*0.5;
  Real drhoy = std::abs(sdatafab(i,j+1,k,URHO) - sdatafab(i,j-1,k,URHO))*0.5;
  Real drhoz = std::abs(sdatafab(i,j,k+1,URHO) - sdatafab(i,j,k-1,URHO))*0.5;
  Real rhop  = sdatafab(i,j,k,URHO) + 1.e-8;
  Real rhofluc = std::sqrt(drhox*drhox + drhoy*drhoy + drhoz*drhoz)/rhop ;

  // rho flux is relative drho in a cell
  //tagfab(i,j,k) = (rhofluc > rhofluc_threshold[level]);

  // always refine close to body (at all levels)
  // use ghost points to refine

  if (nt > 0) {
    if (ibfab(i,j,k,1)) {
      for (int ii = -1; ii <= 1; ii++) {
        for (int jj = -1; jj <= 1; jj++) {
          for (int kk = -1; kk <= 1; kk++) {
          tagfab(i+ii,j+jj,k+kk) = true;
          }
        }
      }
    }
  }
  

}
//////////////////////////// Boundary conditions ///////////////////////////////
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
bcnormal(const amrex::Real x[AMREX_SPACEDIM], amrex::Real dratio, const amrex::Real s_int[ProbClosures::NCONS],
         const amrex::Real s_refl[ProbClosures::NCONS], amrex::Real s_ext[ProbClosures::NCONS],
         const int idir, const int sgn, const amrex::Real time,
         amrex::GeometryData const& /*geomdata*/,  ProbClosures const& closures, ProbParm const& pparm)  
{
  const int UFS  = ProbClosures::UFS;
  const int UMX  = ProbClosures::UMX;
  const int UMY  = ProbClosures::UMY;
  const int UMZ  = ProbClosures::UMZ;
  const int UET  = ProbClosures::UET;
  const int face = (idir+1)*sgn;  

  switch(face)
  {    
    case  -1:  // EAST
      break; 
    case   1:  // WEST
      // inflow      
      s_ext[UMX]  = pparm.rho_oo * pparm.u_oo;
      s_ext[UMY]  = 0.0;
      s_ext[UMZ]  = 0.0;
      s_ext[UET]  = pparm.rho_oo *pparm.eint_oo + pparm.kin_oo;      
      
      s_ext[UFS + N2_ID]  =  0.0;  
      s_ext[UFS + CO2_ID] =  pparm.rho_oo*pparm.Yco2_oo;  
      s_ext[UFS + AR_ID]  =  pparm.rho_oo*pparm.Yar_oo;  
      
      break;  
    default:
      break;
  }


}

///////////////////////// IBM USER CLASS
//// follows ib_walltypes
template < typename param, typename cls_t>
class ibm_user_t
{
  private:
  
  public:

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static void compute_surfIB(const Array1D<Real,0,AMREX_SPACEDIM-1>& xyz,const Array1D<Real,0,AMREX_SPACEDIM-1>& norm,
      const Array1D<Real,0,AMREX_SPACEDIM-1>& t1,const Array1D<Real,0,AMREX_SPACEDIM-1>& t2,
      Array2D<Real,0,ibm_eorder+1,0,cls_t::NPRIM-1>& q, int type_solid_bc, const cls_t* cls) {

        
      // slip velocity (in local coordinates) normal pointing OUT
      q(1,cls_t::QU) = 0.0_rt; // un
      q(1,cls_t::QV) = q(2,cls_t::QV); // ut1
      q(1,cls_t::QW) = q(2,cls_t::QW); // ut2

      // zerograd pressure    
      q(1,cls_t::QPRES) = q(2,cls_t::QPRES); 
      // zerograd temperature
      q(1,cls_t::QT)    = q(2,cls_t::QT);

      // zerograd species (normalise so sum(Y)=1)
#if NUM_SPECIES > 1          
      Real Yw[NUM_SPECIES] ={0.0};
      Real sumY = 0.0;
      for (int n = 0; n < NUM_SPECIES; ++n) {
        Yw[n]   =  q(2,cls_t::QFS+n);
        sumY += Yw[n];
      }
      for (int n = 0; n < NUM_SPECIES; ++n) { 
        q(1,cls_t::QFS+n)   =  Yw[n]/sumY;
      }
#endif                      
     
      // locate the SRP
      
      const Real xjet = xyz(0) - param::xsrp;
      const Real yjet = xyz(1) - param::ysrp;
      const Real zjet = xyz(2) - param::zsrp;
      const Real Rjet=sqrt(yjet*yjet + zjet*zjet);

      bool isjet = (std::fabs(xjet) < 0.01) && (Rjet < param::Rsrp);
      if (isjet)
      {
        q(1,cls_t::QU)    = param::u_srp; 
        q(1,cls_t::QPRES) = param::Pt; 
        q(1,cls_t::QT)    = param::Tt;         
        // SRP-jet composition  (pure Nitrogen)
        Real Yjet[NUM_SPECIES] ={0.0};
        Yjet[N2_ID] = 1.0; 
        for (int n = 0; n < NUM_SPECIES; ++n) {
          q(1,cls_t::QFS+n) = Yjet[n];         
        }
      }

              // std::cout << " Pt = " << param::Pt  << std::endl; 
        // std::cout << " Tt = " << param::Tt << std::endl; 
        // std::cout << " ut = " << param::u_srp << std::endl;         
        // printf(" norm %f %f %f \n", norm(0),norm(1),norm(2));
        // printf(" xyz  %f %f %f \n", xyz(0),xyz(1),xyz(2));
        // std::cout << " P0 = " << param::P0  << std::endl; 
        // std::cout << " T0 = " << param::T0  << std::endl; 
        // std::cout << " T0srp = " << param::T0srp  << std::endl;                 




    }
};    


/////////////////////////////////////////

}
#endif
