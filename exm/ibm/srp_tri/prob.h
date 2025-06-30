#ifndef CNS_PROB_H_
#define CNS_PROB_H_

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_ParmParse.H>

#include <Closures.h>
#include <RHS.h>

#include <eib.h>
#include <Constants.h>
#include <NozzleFunctions.h>
#include <ib_walltypes.h>

#include <numbers>

using namespace amrex;
using namespace universal_constants;

namespace PROB {

// constants  
static constexpr Real Mach     = 4.6;              // bulk Mach number
static constexpr Real Mw       = 28.96e-3;           // Molecular weight [CO2]
static constexpr Real gam      = 1.4;              // Adiabatic coefficient 
static constexpr Real Rgas     = gas_constant/Mw;  // gas constant
static constexpr Real Cv       = Rgas/(gam - 1.0);
static constexpr Real Cp       = gam*Cv;

static constexpr bool srp_on   = false;
static constexpr int ibm_eorder=1;
  
//////////////////////////// Physical modelling ////////////////////////////////
struct ProbParm
{ 

  //------ values obtained from calc_nozzle.py  (SI units)---------
  // based on Mach number 4.6 and gamma 1.4
  // and free-stream conditions P_oo = 574.56 [Pa] and T_oo = 65.0 [K]
  // Staganation conditions
  static constexpr Real P0  =  188222.90041860205 ;
  static constexpr Real T0  =  340.07999999999987 ;
  // Free-stream conditions
  static constexpr Real p_oo  =  574.56 ;
  static constexpr Real T_oo  =  65.0 ;
  static constexpr Real rho_oo  =  0.030799249530956845 ;
  static constexpr Real c_oo  =  161.60754932861272 ;
  static constexpr Real u_oo  =  743.3947269116185 ;
  static constexpr Real eint_oo  =  46637.50000000001 ;
  static constexpr Real kin_oo  =  8510.382719999996 ;
  // SRP Staganation conditions
  static constexpr Real P0srp  =  4437901.4399999995 ;
  static constexpr Real T0srp  =  347.09999999999997 ;
  // SRP throat conditions (choked flow)
  static constexpr Real Pt  =  2344462.5064358213 ;
  static constexpr Real Tt  =  289.25 ;
  static constexpr Real u_srp  =  340.9114987793753 ;
  // SRP nozzle exit pressure (based on Mach number 2.94 )
  static constexpr Real Pe_srp  =  132227.3346346078 ;
 //-----------------------------------------------------------------

  // Wall Temperature
  static constexpr Real Twall = 2.0*T_oo;  // 130 K

  // centre of mass of probe (approx)  
  static constexpr Real x0 = 0.3, y0 = 0.275, z0 = 0.275;

  // number if nozzles and position  (labeled clockwise from top)
  static constexpr int nozzles = 3;  
  static constexpr Real xsrp[nozzles]= { 0.299,0.299,0.299}; 
  static constexpr Real ysrp[nozzles]= { 0.307,0.259,0.259}; 
  static constexpr Real zsrp[nozzles]= { 0.275,0.302,0.248}; 


  // area of nozzle
  static constexpr Real Rsrp= 0.0028; // nozzle radius ( ~ 3 mm)
  static constexpr Real Asrp= std::numbers::pi*Rsrp*Rsrp;  
};

//  parameters for viscous solver and conductivity/viscosity
struct methodparm_t {

  public:

  static constexpr int  order = 2;                  // order numerical scheme   
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
};

// CLOSURES (perfect gas + Sutherland's)
using ProbClosures = closures_dt< indicies_t, transport_suth_t, calorifically_perfect_gas_t<indicies_t> >;


// NUMERICAL SCHEME + EQNS TO SOLVE   (Riemann                  
typedef rhs_dt<riemann_t<false, ProbClosures>, no_diffusive_t, no_source_t > ProbRHS;

//using ProbRHS = rhs_dt< riemann_t<false, ProbClosures>, viscous_t<methodparm_t, ProbClosures>, reactor_t<ProbClosures> >;

// declaration of user-specific ibm (see bottom of the file for definition)
template < typename param, typename cls_t > class ibm_user_t; 

// IBM templates
typedef ibm_user_t<ProbParm,ProbClosures> TypeWall;
typedef eib_t<TypeWall,ibmparm_t,ProbClosures> ProbIB;


void inline inputs() {
  
  amrex::Print() << " ****** Starting ... ******* " <<  std::endl;
  amrex::Print() << " SRP tri-nozzle  " <<  std::endl;
  amrex::Print() << " adapated from Jiaye 2022 " <<  std::endl;  
  amrex::Print() << " **************************** " <<  std::endl;

  ProbParm  data;

  amrex::Print() << " (ATM) P= " << data.p_oo << " T=" << data.T_oo << std::endl;
  amrex::Print() << " (ATM) P0= " << data.P0 << " T0=" << data.T0 << std::endl;
  amrex::Print() << " (ATM) u= " << data.u_oo << std::endl;

  amrex::Print() << " (INLET SRP) P= " << data.Pt << " T=" << data.Tt << std::endl;
  amrex::Print() << " (INLET SRP) P0= " << data.P0srp << " T0=" << data.T0srp << std::endl;
  amrex::Print() << " (INLET SRP) u= " << data.u_srp << std::endl;
  amrex::Print() << " (expected OUTLET SRP) Pe= " << data.Pe_srp << std::endl;
  amrex::Print() << " Pe/P0 " << data.Pe_srp/data.P0 << std::endl;
  amrex::Print() << " **************************** " <<  std::endl;


}

//////////////////////////// Initial conditions ////////////////////////////////
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void prob_initdata (int i, int j, int k, amrex::Array4<amrex::Real> const& state,
      amrex::GeometryData const& geomdata, ProbClosures const& cls, ProbParm const& pparm) {
  
  const Real* prob_lo = geomdata.ProbLo();  
  const Real* dx      = geomdata.CellSize();

  Real x = prob_lo[0] + (i+0.5_rt)*dx[0];
  // Real y = prob_lo[1] + (j+0.5_rt)*dx[1];
  // Real z = prob_lo[2] + (k+0.5_rt)*dx[2];

  // local vars
  Real rhot,eint,u[3]={0.0};

  // initial state
  rhot =  pparm.rho_oo;
  u[0] =  pparm.u_oo;
  if (x > 0.25) u[0] = 0.0;
  eint =  pparm.rho_oo*pparm.eint_oo;

  state(i, j, k, cls.URHO) = rhot ;
  state(i, j, k, cls.UMX)  = rhot * u[0];
  state(i, j, k, cls.UMY)  = rhot * u[1];
  state(i, j, k, cls.UMZ)  = rhot * u[2];
  state(i, j, k, cls.UET)  = eint + Real(0.5) * rhot * u[0] * u[0] ; 
  
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE 
void user_tagging(int i, int j, int k, int nt, auto& tagfab, const auto &sdatafab, 
                  const Array4<bool>&ibfab, const auto& geomdata, 
                  const ProbParm& pparm , int level) {

  const Real* dx  = geomdata.CellSize();
  const Real* prob_lo = geomdata.ProbLo();  
  const Real x = prob_lo[0] + (i+0.5_rt)*dx[0];
  const Real y = prob_lo[1] + (j+0.5_rt)*dx[1];
  const Real z = prob_lo[2] + (k+0.5_rt)*dx[2];

  bool refine = false;

  // close to domain Box [ 0.25-0.35,0.20-0.35,0.20-0.35]
  // RULE 1: refine within box (only levels 0 and 1)

  //const bool refinebox =  
  //((x > 0.25) && (x < 0.35 ) && (y > 0.15) && (y < 0.4 ) && (z > 0.15) && (z < 0.4 ));
  //if (level < 2) refine = refinebox;

  // RULE 2: refine if gradients of density large
  // initialize thresholds at all levels
  const Real rhofluc_threshold[6] = {0.02_rt,0.05_rt,0.1_rt,0.1_rt,0.1_rt,0.7_rt};
  int URHO = ProbClosures::URHO; 
  // refine close to grads of rho
  Real drhox = std::abs(sdatafab(i+1,j,k,URHO) - sdatafab(i-1,j,k,URHO))*0.5;
  Real drhoy = std::abs(sdatafab(i,j+1,k,URHO) - sdatafab(i,j-1,k,URHO))*0.5;
  Real drhoz = std::abs(sdatafab(i,j,k+1,URHO) - sdatafab(i,j,k-1,URHO))*0.5;
  Real rhop  = sdatafab(i,j,k,URHO);
  Real rhofluc = std::sqrt(drhox*drhox + drhoy*drhoy + drhoz*drhoz)/rhop ;  
  refine = (rhofluc > rhofluc_threshold[level]) || refine;

  // RULE 3: refine close to body and within box) 
  //  (only after nt> 0)
  // if ((nt > 0) && refinebox) {
  //   bool foundGP= false;
  //   for (int ii = -1; ii <= 1; ii++) {
  //     for (int jj = -1; jj <= 1; jj++) {
  //       for (int kk = -1; kk <= 1; kk++) {          
  //         foundGP =   ibfab(i,j,k,1) || foundGP;               
  //       }
  //     }
  //   }
  //   refine = foundGP || refine;
  // }

  // RULE 4: refine close to nozzles
  for (int nozz=0;nozz<pparm.nozzles;nozz++) {
    Real xjet = x - pparm.xsrp[nozz];
    Real yjet = y - pparm.ysrp[nozz];
    Real zjet = z - pparm.zsrp[nozz];
    // distance to centre nozzle
    Real dis = sqrt(xjet*xjet + yjet*yjet + zjet*zjet);
    refine = (dis < 0.01) || refine;  
  }

  tagfab(i,j,k) = refine;
  

}
//////////////////////////// Boundary conditions ///////////////////////////////
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
bcnormal(const amrex::Real x[AMREX_SPACEDIM], amrex::Real dratio, const amrex::Real s_int[ProbClosures::NCONS],
         const amrex::Real s_refl[ProbClosures::NCONS], amrex::Real s_ext[ProbClosures::NCONS],
         const int idir, const int sgn, const amrex::Real time,
         amrex::GeometryData const& /*geomdata*/,  ProbClosures const& closures, ProbParm const& pparm)  
{
  const int URHO  = ProbClosures::URHO;
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
      s_ext[URHO] = pparm.rho_oo;            
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
      Array2D<Real,0,ibm_eorder+1,0,cls_t::NPRIM-1>& q, const cls_t* cls) {

        
      // slip velocity (in local coordinates) normal pointing OUT
      q(1,cls_t::QU) = 0.0_rt; // un
      q(1,cls_t::QV) = q(2,cls_t::QV); // ut1
      q(1,cls_t::QW) = q(2,cls_t::QW); // ut2

      // zerograd pressure    
      q(1,cls_t::QPRES) = q(2,cls_t::QPRES); 
      // zerograd temperature
      q(1,cls_t::QT)    = q(2,cls_t::QT);

      // activate SRP 
      if (srp_on)
      {
        // locate the SRP ------- (looking in the three nozzles)
        // NOTE: all xyz points are in the surface
        for (int nozz=0;nozz<param::nozzles;nozz++) 
        {
          // vector to centre of nozzle  xjet
          Real xjet = xyz(0) - param::xsrp[nozz];
          Real yjet = xyz(1) - param::ysrp[nozz];
          Real zjet = xyz(2) - param::zsrp[nozz];
          // distance to centre nozlle
          Real Rjet = sqrt(xjet*xjet + yjet*yjet + zjet*zjet);
          // detects points close to centre nozzle        
          if (Rjet < param::Rsrp)
          {
            q(1,cls_t::QU) = param::u_srp; 
            q(1,cls_t::QV) = 0.0; 
            q(1,cls_t::QW) = 0.0; 
            q(1,cls_t::QPRES) = param::Pt; 
            q(1,cls_t::QT)    = param::Tt;         
          }
      }
      }
      // -----------

    }
};    


/////////////////////////////////////////

}
#endif
