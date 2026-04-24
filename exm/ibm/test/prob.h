#ifndef CNS_PROB_H_
#define CNS_PROB_H_

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_ParmParse.H>

#include <Closures.h>
#include <RHS.h>

#include <ibm_solver.h>
#include <Constants.h>
#include <NozzleFunctions_jiaye.h>
#include <ibm_walltypes.h>
#include <numbers>

using namespace amrex;
using namespace universal_constants;

namespace PROB {

// constants  
static constexpr Real Mw       = 28.96e-3;         // Molecular weight [air]
static constexpr Real gam      = 1.4;              // Adiabatic coefficient 
static constexpr Real Rgas     = gas_constant/Mw;  // gas constant
static constexpr Real Cv       = Rgas/(gam - 1.0);
static constexpr Real Cp       = gam*Cv;

static constexpr bool srp_on   = true;
static constexpr int ibm_eorder= 1;
  
//////////////////////////// Physical modelling ////////////////////////////////
struct ProbParm
{ 

  // freee-stream conditions  (NASA tunnel)
  static constexpr Real Ma_oo   = 4.6;                    // [1]  free-stream Mach number
  static constexpr Real p_oo    = 574.56;                 // [Pa] free-stream Pressure   
  static constexpr Real T_oo    = 65.0;                   // [K]  free-stream Temperature 

  // Wall Temperature
  static constexpr Real Twall   = 4.0*T_oo;
  
  // free-stream stagnation pressure and temperature
  //static constexpr Real P0 = nozzle_functions::Pstag(p_oo,Ma_oo,gam);
  //static constexpr Real T0 = nozzle_functions::Tstag(T_oo,Ma_oo,gam);

  // free-stream stagnation pressure and temperature (from calc_nozzle.py)
  // these lines replace the above values in clang compilers
  static constexpr Real P0  =  188222.90041860205 ;
  static constexpr Real T0  =  340.07999999999987 ;

  // free-stream conditions (from calc_nozzle.py) 
  static constexpr Real rho_oo  =  0.030799249530956845 ;
  static constexpr Real c_oo  =  161.60754932861272 ;
  static constexpr Real u_oo  =  743.3947269116185 ;
  static constexpr Real eint_oo  =  46637.50000000001 ;
  static constexpr Real kin_oo  =  8510.382719999996 ; 


  // geo paramaters
  // model attitude 
  static constexpr Real x0             = 0.0;                 // center On the base of FOREbody CONE/m
  static constexpr Real y0             = 0.0;                 // center On the base of FOREbody CONE/m
  static constexpr Real z0             = 0.0;                 // center On the base of FOREbody CONE/m
  static constexpr Real Geo_AoA        = 0.0;                 // angle of attack/deg
  //static constexpr Real Geo_RollA      = 0.0;                 // roll angle/deg
  //static constexpr Real Geo_YawA       = 0.0;                 // yaw angle/deg
  // model paramaters
  static constexpr Real Diameter_Cone  = 0.127;                // forebody cone diameter/m
  static constexpr Real Radius_Cone    = 0.5*Diameter_Cone;    // forebody cone radius/m
  static constexpr Real Theta_Cone_Deg = 70.0;                 // forebody cone apex angle/deg
  static constexpr Real Len_Aftbody    = 0.24257;              // aftbody cylinder length/m
  static constexpr Real Nz_Exit_To_Geo = 0.1;                  // ratio of nozzle exit diameter to forebody cone diameter
  // static constexpr Real Height_Cone    = Radius_Cone/std::tan(Theta_Cone_Deg*std::numbers::pi/180);  // forebody cone height/m

  static constexpr Real Height_Cone  =  0.02311210987590386 ;

  // Nozzle paramaters
  // number of nozzles
  static constexpr int  Nozzle_Count   = 4;
  // nozzle conditions
  static constexpr Real Ma_srp         = 2.94018;                                         // nozzle exit Mach number
  static constexpr Real P0srp          = 7724.0*p_oo;                                     // nozzle stagnation pressure
  static constexpr Real T0srp          = 5.34*T_oo;                                       // nozzle stagnation tempurture
  static constexpr Real R_Nozzle_Exit  = Nz_Exit_To_Geo*Radius_Cone;                      // radius of nozzle exit/m
  //static constexpr Real R_Plenum       = 0.0;                                             // radius of nozzle plenum/m
  //static constexpr Real Divergence_Deg = 15.0;                                            // nozzle expansion angle/deg
  //static constexpr Real Convergence_Deg= 15.0;                                            // nozzle Compression angle/deg

  static constexpr int Nozzle_Switch[Nozzle_Count] = {1,0,0,0};  // 1 using this nozzle; 0 block this nozzle
  static constexpr int Nozzle_Type[Nozzle_Count]   = {1,0,0,0};  // 1 apply bc at thorat; 0 apply bc at exit; 2 apply bc at stagnation condition

   // location of Centre of nozzle exit, for x0, y0, z0 = 0, aoa = 0
  static constexpr std::array<std::array<Real, AMREX_SPACEDIM>, Nozzle_Count> Nozzle_Exit_Centre = {
    // Central Nozzle location
    std::array<Real, AMREX_SPACEDIM>({   -0.0208 ,       0.0 ,     0.0 }),      //[ 0    , 0/1 ]
    // Tri Nozzle at Half radius
    std::array<Real, AMREX_SPACEDIM>({ -0.009245 ,   0.03175 ,     0.0 }),      //[ 0    , 1/2 ]   
    std::array<Real, AMREX_SPACEDIM>({ -0.009245 , -0.015875 ,  0.0275 }),      //[ 120  , 1/2 ]
    std::array<Real, AMREX_SPACEDIM>({ -0.009245 , -0.015875 , -0.0275 }) };    //[ 240  , 1/2 ]
    // Tri Nozzle at 1/4 radius
    //std::array<Real, AMREX_SPACEDIM>({ −0.01502 ,   0.015875 ,     0.0 }),      //[ 0    , 1/4 ]
    //std::array<Real, AMREX_SPACEDIM>({ −0.01502 , −0.0079375 ,  0.01375 }),     //[ 120  , 1/4 ]
    //std::array<Real, AMREX_SPACEDIM>({ −0.01502 , −0.0079375 , -0.01375 }) };   //[ 240  , 1/4 ]
    // Tri Nozzle at 3/4 radius
    //std::array<Real, AMREX_SPACEDIM>({ −0.003467 ,   0.047625 ,     0.0 }),     //[ 0    , 3/4 ]
    //std::array<Real, AMREX_SPACEDIM>({ −0.003467 , −0.0238125 ,  0.04125 }),    //[ 120  , 3/4 ]
    //std::array<Real, AMREX_SPACEDIM>({ −0.003467 , −0.0238125 , -0.04125 }) };  //[ 240  , 3/4 ]
    // Quad Nozzle at Half radius
    //std::array<Real, AMREX_SPACEDIM>({ -0.009245 ,   0.03175 ,      0.0 }),     //[ 0    , 1/2 ]
    //std::array<Real, AMREX_SPACEDIM>({ -0.009245 ,   0.0     ,  0.03175 }),     //[ 90   , 1/2 ]
    //std::array<Real, AMREX_SPACEDIM>({ -0.009245 ,  -0.03175 ,      0.0 }),     //[ 180  , 1/2 ]
    //std::array<Real, AMREX_SPACEDIM>({ -0.009245 ,   0.0     , -0.03175 }) };   //[ 270  , 1/2 ]
    // Hexa Nozzle at Half radius
    //std::array<Real, AMREX_SPACEDIM>({ -0.009245 ,   0.03175 ,     0.0 }),      //[ 0    , 1/2 ]
    //std::array<Real, AMREX_SPACEDIM>({ -0.009245 ,  0.015875 ,  0.0275 }),      //[ 60   , 1/2 ]
    //std::array<Real, AMREX_SPACEDIM>({ -0.009245 , -0.015875 ,  0.0275 }),      //[ 120  , 1/2 ]
    //std::array<Real, AMREX_SPACEDIM>({ -0.009245 ,   0.03175 ,     0.0 }),      //[ 180  , 1/2 ]
    //std::array<Real, AMREX_SPACEDIM>({ -0.009245 , -0.015875 , -0.0275 }) };    //[ 240  , 1/2 ]
    //std::array<Real, AMREX_SPACEDIM>({ -0.009245 ,  0.015875 , -0.0275 }) };    //[ 300  , 1/2 ]

  // nozzle orientation for aoa = 0, roll angle = 0
  static constexpr std::array<std::array<Real, AMREX_SPACEDIM>, Nozzle_Count> Nozzle_Orientation = {
    std::array<Real, AMREX_SPACEDIM>({     -1.0 ,      0.0 ,     0.0 }),
    // Nozzle orientation parallel to the Geo axis
    std::array<Real, AMREX_SPACEDIM>({     -1.0 ,      0.0 ,     0.0 }),
    std::array<Real, AMREX_SPACEDIM>({     -1.0 ,      0.0 ,     0.0 }),
    std::array<Real, AMREX_SPACEDIM>({     -1.0 ,      0.0 ,     0.0 }) };
    // Nozzle orientation perpendicular to surface
    //std::array<Real, AMREX_SPACEDIM>({-std::sin(Theta_Cone_Deg*std::numbers::pi/180) ,                                    
    //                                   std::cos(Theta_Cone_Deg*std::numbers::pi/180) , 0.0 }),
    //std::array<Real, AMREX_SPACEDIM>({-std::sin(Theta_Cone_Deg*std::numbers::pi/180) ,                                    
    //                                   std::cos(Theta_Cone_Deg*std::numbers::pi/180) * std::cos(120*std::numbers::pi/180), 
    //                                   std::cos(Theta_Cone_Deg*std::numbers::pi/180) * std::sin(120*std::numbers::pi/180) }),
    //std::array<Real, AMREX_SPACEDIM>({-std::cos(20*std::numbers::pi/180) ,                                    
    //                                   std::cos(Theta_Cone_Deg*std::numbers::pi/180) * std::cos(240*std::numbers::pi/180), 
    //                                   std::cos(Theta_Cone_Deg*std::numbers::pi/180) * std::sin(240*std::numbers::pi/180) }) }; 

  static constexpr int Nozzles_ON = std::accumulate(Nozzle_Switch, Nozzle_Switch + Nozzle_Count, 0);

  static constexpr nozzle_functions::Nozzle Nozzles[Nozzle_Count] = {
    nozzle_functions::Nozzle(Nozzle_Type[0], Nozzle_Switch[0], Nozzle_Exit_Centre[0], Nozzle_Orientation[0], Ma_srp, P0srp/Nozzles_ON, T0srp, R_Nozzle_Exit),
    nozzle_functions::Nozzle(Nozzle_Type[1], Nozzle_Switch[1], Nozzle_Exit_Centre[1], Nozzle_Orientation[1], Ma_srp, P0srp/Nozzles_ON, T0srp, R_Nozzle_Exit),
    nozzle_functions::Nozzle(Nozzle_Type[2], Nozzle_Switch[2], Nozzle_Exit_Centre[2], Nozzle_Orientation[2], Ma_srp, P0srp/Nozzles_ON, T0srp, R_Nozzle_Exit),
    nozzle_functions::Nozzle(Nozzle_Type[3], Nozzle_Switch[3], Nozzle_Exit_Centre[3], Nozzle_Orientation[3], Ma_srp, P0srp/Nozzles_ON, T0srp, R_Nozzle_Exit) };

};


//  parameters for viscous solver and conductivity/viscosity
struct methodparm_t {

  public:

  static constexpr int  order = 2;                  // order numerical scheme  
  static constexpr bool use_LES = false;
};


// parameters for skew-symmetric method
struct skewparm_t {

  public:

  static constexpr bool dissipation = true;         // no dissipation
  static constexpr int  order       = 4;            // order numerical scheme   
  static constexpr Real C2skew=1.5,C4skew=0.016;    // Skew symmetric values 
};

struct ibmparm_t {

  public:

  static constexpr int  interp_order = 1;
  static constexpr int  extrap_order = ibm_eorder;
  static constexpr Real alpha        = 0.6;      
};

// CLOSURES (perfect gas + Sutherland's)
using ProbClosures = closures_dt< indicies_t, transport_suth_t, calorifically_perfect_gas_t<indicies_t> >;


// NUMERICAL SCHEME + EQNS TO SOLVE   (Riemann                  
typedef rhs_dt<riemann_t<false, ProbClosures>, viscous_t<methodparm_t, ProbClosures>, no_source_t > ProbRHS;

//using ProbRHS = rhs_dt< riemann_t<false, ProbClosures>, viscous_t<methodparm_t, ProbClosures>, reactor_t<ProbClosures> >;

// declaration of user-specific ibm (see bottom of the file for definition)
template < typename param, typename cls_t > class ibm_user_t; 

// IBM templates
typedef ibm_user_t<ProbParm,ProbClosures> TypeWall;
typedef ibm_solver_t<TypeWall,ibmparm_t,ProbClosures> ProbIB;


void inline inputs() {
  
  amrex::Print() << " ****** Starting ... ******* " <<  std::endl;
  amrex::Print() << " SRP tri-nozzle  " <<  std::endl;
  amrex::Print() << " adapated from Jiaye 2022 " <<  std::endl;  
  amrex::Print() << " **************************** " <<  std::endl;

  ProbParm  data;

  amrex::Print() << " (ATM) P = " << data.p_oo << " T=" << data.T_oo << std::endl;
  amrex::Print() << " (ATM) P0= " << data.P0 << " T0=" << data.T0 << std::endl;
  amrex::Print() << " (ATM) u = " << data.u_oo << std::endl;

  amrex::Print() << " (INLET SRP) P_Throat= " << data.Nozzles[0].P_Throat() << " T_Throat=" << data.Nozzles[0].T_Throat() << std::endl;
  amrex::Print() << " (INLET SRP) P_Stag= " << data.Nozzles[0].P_Stag << " T_Stag=" << data.Nozzles[0].T_Stag << std::endl;
  amrex::Print() << " (INLET SRP) V_Throat= " << data.Nozzles[0].V_Throat() << std::endl;
  amrex::Print() << " (expected OUTLET SRP) Mach_Exit= " << data.Nozzles[0].Mach_Exit   << std::endl;
  amrex::Print() << " (expected OUTLET SRP) P_Exit= " << data.Nozzles[0].P_Exit() << std::endl;
  amrex::Print() << " (expected OUTLET SRP) T_Exit= " << data.Nozzles[0].T_Exit() << std::endl;
  amrex::Print() << " Pe/P0= " << data.Nozzles[0].P_Exit()/data.P0 << std::endl;
  amrex::Print() << " ***************************************** " <<  std::endl;


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
  if (x > (pparm.x0 - pparm.Radius_Cone)) {u[0] = 0.0;rhot = 0.5*rhot;}
  eint =  pparm.rho_oo*pparm.eint_oo;

  state(i, j, k, cls.URHO) = rhot;
  state(i, j, k, cls.UMX)  = rhot * u[0];
  state(i, j, k, cls.UMY)  = rhot * u[1];
  state(i, j, k, cls.UMZ)  = rhot * u[2];
  state(i, j, k, cls.UET)  = eint + Real(0.5) * rhot * u[0] * u[0] ; 

}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE 
void user_tagging(int i, int j, int k, int nt, auto& tagfab, const auto &sdatafab, 
                  const Array4<bool>&ibfab, const auto& geomdata, 
                  const ProbParm& pparm , int level) {

  const Real* dx      = geomdata.CellSize();
  const Real* prob_lo = geomdata.ProbLo();  
  const Real x = prob_lo[0] + (i+0.5_rt)*dx[0];
  const Real y = prob_lo[1] + (j+0.5_rt)*dx[1];
  const Real z = prob_lo[2] + (k+0.5_rt)*dx[2];

  bool refine = false;

if (level < 3) {
  // RULE 2: refine if gradients of density large
  // initialize thresholds at all levels
  const Real rhofluc_threshold[6] = {0.1_rt,0.1_rt,0.2_rt,0.2_rt,0.5_rt,0.7_rt};
  int URHO = ProbClosures::URHO; 
  // refine close to grads of rho
  Real drhox = std::abs(sdatafab(i+1,j,k,URHO) - sdatafab(i-1,j,k,URHO))*0.5;
  Real drhoy = std::abs(sdatafab(i,j+1,k,URHO) - sdatafab(i,j-1,k,URHO))*0.5;
  Real drhoz = std::abs(sdatafab(i,j,k+1,URHO) - sdatafab(i,j,k-1,URHO))*0.5;
  Real rhop  = sdatafab(i,j,k,URHO);
  Real rhofluc = std::sqrt(drhox*drhox + drhoy*drhoy + drhoz*drhoz)/rhop ;  
  refine = (rhofluc > rhofluc_threshold[level]) || refine;
}


if (level < 4) {
  // RULE 3: refine close to body //(only after nt> 0)
  if (nt > 0) {
     bool foundGP= false;
     for (int ii = -1; ii <= 1; ii++) {
       for (int jj = -1; jj <= 1; jj++) {
         for (int kk = -1; kk <= 1; kk++) {          
           foundGP =   ibfab(i+ii,j+jj,k+kk,1) || foundGP;               
         }
       }
     }
     refine = foundGP || refine;
  }
}


if (level < 4)
{
  // RULE 4: refine close to Nozzle_Count
  for (int nozz=0; nozz<pparm.Nozzle_Count; nozz++) 
  {
    if (pparm.Nozzles[nozz].Nozzle_Switch == 0) continue;

    Real theta      = pparm.Geo_AoA*std::numbers::pi/180.0;
    Real Cen_jet[3] = {
       (pparm.Nozzles[nozz].Nozzle_Centre()[0])*std::cos(theta) + pparm.Nozzles[nozz].Nozzle_Centre()[1]*std::sin(theta),
      -(pparm.Nozzles[nozz].Nozzle_Centre()[0])*std::sin(theta) + pparm.Nozzles[nozz].Nozzle_Centre()[1]*std::cos(theta),
        pparm.Nozzles[nozz].Nozzle_Centre()[2]};

    Real xjet = x - (Cen_jet[0] + pparm.x0);
    Real yjet = y - (Cen_jet[1] + pparm.y0);
    Real zjet = z - (Cen_jet[2] + pparm.z0);

    //Real cross_x = pparm.Nozzles[nozz].Orientation[1] * zjet - pparm.Nozzles[nozz].Orientation[2] * yjet;
    //Real cross_y = pparm.Nozzles[nozz].Orientation[2] * xjet - pparm.Nozzles[nozz].Orientation[0] * zjet;;
    //Real cross_z = pparm.Nozzles[nozz].Orientation[0] * yjet - pparm.Nozzles[nozz].Orientation[1] * xjet;

    //Real cross_norm = std::sqrt(cross_x*cross_x + cross_y*cross_y + cross_z*cross_z);
    //Real Distance_to_NozzleAxis = std::abs(cross_norm/pparm.Nozzles[nozz].Orientation_Norm());
    //refine = (Distance_to_NozzleAxis < 1.2*pparm.Nozzles[nozz].R_Exit) || refine; 

    // distance to centre nozzle volume
    Real dis = xjet*xjet + yjet*yjet + zjet*zjet;
    const Real containment = 1.1;

    refine = (dis < containment*pparm.Nozzles[nozz].Nozzle_Vol()) || refine;     
  }
}


if (level < 4) 
{
  // RULE 5: refine within cylinder enclosing the forebody cone
  Real theta      = pparm.Geo_AoA*std::numbers::pi/180.0;
  Real x_resolved = (x - pparm.x0)*std::cos(theta) - (y - pparm.y0)*std::sin(theta);
  Real y_resolved = (x - pparm.x0)*std::sin(theta) + (y - pparm.y0)*std::cos(theta);
  Real z_resolved = (z - pparm.z0);
  Real r_resloved = std::sqrt(y_resolved*y_resolved + z_resolved*z_resolved);

  bool refinebox = ( 
    ( x_resolved > -1.5*pparm.Height_Cone) && (x_resolved < 0.5*pparm.Height_Cone) && 
      r_resloved <  1.1*pparm.Radius_Cone);                                                   //enclosing the forebody cone by cylinder
    //(y_resolved > -1.1*pparm.Radius_Cone ) && (y_resolved < 1.1*pparm.Radius_Cone ) &&  
    //(z_resolved > -1.1*pparm.Radius_Cone ) && (z_resolved < 1.1*pparm.Radius_Cone ) );      //enclosing the forebody cone by Rectangular box
  if (refinebox) {
  const Real proximity = 0.95;
  refinebox = ((x_resolved + proximity*pparm.Height_Cone)/r_resloved) <= std::tan((90-pparm.Theta_Cone_Deg)*std::numbers::pi/180.0); } 

  refine = refinebox || refine;
}

  //if (x > 0.2) refine=false;
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
      q(1,cls_t::QV) = 0.0_rt; //q(2,cls_t::QV); // ut1
      q(1,cls_t::QW) = 0.0_rt; //q(2,cls_t::QW); // ut2

      // zerograd pressure    
      q(1,cls_t::QPRES) = q(2,cls_t::QPRES); 
      // zerograd temperature
      q(1,cls_t::QT)    = param::Twall;

      // activate SRP 
      if (srp_on)
      {
        // locate the SRP ------- (loop in Nozzle_Count)
        // NOTE: all xyz points are in the surface
        for (int nozz=0; nozz < param::Nozzle_Count; ++nozz) 
        {
          if ( !param::Nozzles[nozz].Nozzle_Switch ) continue;

          Real theta = param::Geo_AoA*std::numbers::pi/180.0;
          Real Cen_jet[3] = {
              (param::Nozzles[nozz].BC_Applied_Centre()[0])*std::cos(theta) + param::Nozzles[nozz].BC_Applied_Centre()[1]*std::sin(theta),
             -(param::Nozzles[nozz].BC_Applied_Centre()[0])*std::sin(theta) + param::Nozzles[nozz].BC_Applied_Centre()[1]*std::cos(theta),
               param::Nozzles[nozz].BC_Applied_Centre()[2]};

          Real xjet = xyz(0) - (Cen_jet[0] + param::x0);
          Real yjet = xyz(1) - (Cen_jet[1] + param::y0);
          Real zjet = xyz(2) - (Cen_jet[2] + param::z0);
          // distance to centre nozzle
          Real Rjet = std::sqrt(xjet*xjet + yjet*yjet + zjet*zjet);
          // detects points close to centre nozzle
          if (Rjet < param::Nozzles[nozz].R_BC_Applied())
          {
            q(1,cls_t::QU)    = param::Nozzles[nozz].V_BC_Applied(); 
            q(1,cls_t::QV)    = 0.0;
            q(1,cls_t::QW)    = 0.0; 
            q(1,cls_t::QPRES) = param::Nozzles[nozz].P_BC_Applied(); 
            q(1,cls_t::QT)    = param::Nozzles[nozz].T_BC_Applied();         
          }
      }
      }

    }
};    


/////////////////////////////////////////
}
#endif
