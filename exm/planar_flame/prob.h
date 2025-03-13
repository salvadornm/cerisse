#ifndef CNS_PROB_H_
#define CNS_PROB_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_ParmParse.H>
#include <PelePhysics.H>
#include <ReactorBase.H>

#include "Constants.h"
#include "CombustionFunctions.h"

#include "Closures.h"
#include "RHS.h"


namespace PROB {


static constexpr Real Phi   = 0.5;    // Equivalence Ratio
static combustion_functions::Moles_Mixture mix = combustion_functions::compute_MixMol_fromEquivalenceRatio_H2Air(Phi);    // mixture from equivalence Ratio


// SPECIES  (Li and Dryer mechanism)
// H2 O2 H2O H O OH HO2 H2O2 N2
// Fresh gases
//
//       temperature   298 K
//          pressure   1.0132e+05 Pa
//           density   0.98933 kg/m^3
//  mean mol. weight   24.192 kg/kmol
//
//                          1 kg             1 kmol     
//                     ---------------   ---------------
//          enthalpy           -130.24           -3150.7  J
//   internal energy       -1.0255e+05       -2.4809e+06  J
//           entropy            7880.3        1.9064e+05  J/K
//
//                      mass frac. Y      mole frac. X     chem. pot. / RT
//                     ---------------   ---------------   ---------------
//                H2          0.014468           0.17361           -16.827
//                O2           0.22963           0.17361           -25.822
//                N2            0.7559           0.65278           -22.838


// Burning composition obtained with Cantera equivalence ratio 0.5
// using  cerisse/tools/combustion/equil_from_equivalenceratio.py
// 
//       temperature   1644.8 K
//          pressure   1.0132e+05 Pa
//           density   0.19626 kg/m^3
//  mean mol. weight   26.489 kg/kmol
//
//                         1 kg             1 kmol     
//                     ---------------   ---------------
//          enthalpy           -130.24           -3449.9  J
//   internal energy       -5.1641e+05       -1.3679e+07  J
//           entropy            9624.9        2.5496e+05  J/K
//
//                      mass frac. Y      mole frac. X     chem. pot. / RT
//                     ---------------   ---------------   ---------------
//                H2        5.1557e-07        6.7744e-06           -30.764
//                O2           0.11471          0.094966           -30.432
//               H2O           0.12917           0.18993            -45.98
//                 H        8.5022e-09        2.2343e-07           -15.382
//                 O        4.5034e-06        7.4562e-06           -15.216
//                OH        0.00021046        0.00032781           -30.598
//               HO2        4.7965e-07        3.8495e-07           -45.814
//              H2O2        3.3189e-08        2.5847e-08           -61.196
//                N2            0.7559           0.71477           -26.626


// problem parameters
struct ProbParm {

  // unburn gases
  Real rho_u = 0.98933;               // density  [kg/m^3]
  Real T_u   = 298;                   // temperature [K] 
  Real p_u   = 1.0132e+05;            // pressure [Pa]  (5 atm)
  Real e_u   = -1.0255e+05;           // internal energy [J/kg]  
 
  combustion_functions::Y_Mixture massfr = combustion_functions::mol2Y_H2Air(mix);  // unburn mass fractions
  //GpuArray<Real, NUM_SPECIES> Y_u = { massfr.H2, massfr.O2, 0.,0., 0., 0., 0.,0., massfr.N2};  // mass fractions [-] 
  GpuArray<Real, NUM_SPECIES> Y_u = { 0.014468, 0.22963 , 0.,0., 0., 0., 0.,0., 0.7559};  // mass fractions [-] 

  // burn gases
  Real rho_b = 0.19626;                 // density  [kg/m^3]
  Real T_b   = 1644.8;                  // temperature [K]  
  Real p_b   = p_u;                     // pressure [Pa]  (1 atm)   
  Real e_b   = -5.1641e+05 ;            // internal energy [J/kg]
  Real drho = rho_b - rho_u;

 
  // butn mass fractions
  GpuArray<Real, NUM_SPECIES> Y_b = {5.1557e-07,0.11471,0.12917,8.5022e-09,4.5034e-06,0.00021046,4.7965e-07,3.3189e-08,0.7559};

  // geometrical parameters                                     
  Real Lx     =   0.04;  // half-width domain
  Real Ly     =   0.06;
  Real Yflame =   0.5*Ly;

  Real SL     =  0.49; // estimated burning velocity
  Real lf     =  423e-6; // estimated flame thickness  (423 microns)
  // velocity
  Real u_u   = 1.4*SL;                        // velocity [m/s]
  Real u_b   = (SL*drho + rho_u*u_u)/rho_b;   // velocity [m/s]  (RH)
  
  Real mflow =  rho_u*u_u;  // flow rate (per area)

};

// numerical method parameters
struct methodparm_t {

  public:

  static constexpr int  order = 2;        

};



using ProbClosures = closures_dt< indicies_t, transport_Pele_t, multispecies_pele_gas_t<indicies_t> >;

using ProbRHS = rhs_dt< riemann_t<false, ProbClosures>, viscous_t<methodparm_t, ProbClosures>, reactor_t<ProbClosures> >;
//using ProbRHS = rhs_dt< riemann_t<false, ProbClosures>, viscous_t<methodparm_t, ProbClosures>, no_source_t>;
//using ProbRHS = rhs_dt< riemann_t<false, ProbClosures>, no_diffusive_t, reactor_t<ProbClosures> >;



void inline inputs() {
  // ParmParse pp;  
  amrex::Print() << " ****** Starting  *******" <<  std::endl;
  amrex::Print() << " Planar H2-Air Flame (February 2025)  " <<  std::endl;
  amrex::Print() << " Nspecies = " << NUM_SPECIES << std::endl;
  amrex::Print() << " H2 O2 H2O H O OH HO2 H2O2 N2 " << std::endl;
  amrex::Print() << " ******           *******" << std::endl;
}

// initial condition
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void prob_initdata(
    int i, int j, int k, Array4<Real> const &state,
    GeometryData const &geomdata, ProbClosures const &cls,
    ProbParm const &prob_parm) {
  const Real *prob_lo = geomdata.ProbLo();
  const Real *prob_hi = geomdata.ProbHi();
  const Real *dx = geomdata.CellSize();

  Real x = prob_lo[0] + (i + Real(0.5)) * dx[0];
  Real y = prob_lo[1] + (j + Real(0.5)) * dx[1];
  
  Real rhot, vxt,et,Tt,Pt;
  Real Yt[NUM_SPECIES];const Real *Yu;const Real *Yb;

  Real wavelength = prob_parm.Lx;
  Real wave= 2.0*std::numbers::pi/wavelength;
  // pertubation and frequency 
  const Real ypertur = 2*prob_parm.lf; // same of flame thickness (orig 0.04)
  const Real Nwaves = 1;

  Real yinterf =  prob_parm.Yflame + ypertur*cos(Nwaves*wave*x);
  // for (int n=1; n <= Nwaves; ++n){
  //   yinterf += ypertur/n*cos(Nwaves*wave*x*n);
  // }
  
  // smooth profile over lf 
  const Real smooth=2; // the smaller the wider the spread (smooth 5, spread over -lf:lf)
  Real H = (y-yinterf)/prob_parm.lf;
  Real fburn = 0.5*tanh(smooth*H) + 0.5; // fburn   0:unburn  1:burn

  fburn= 0; //switch off flame (for spark)

  Yu  = prob_parm.Y_u.data();
  Yb  = prob_parm.Y_b.data();
  // compute 
  vxt  = prob_parm.u_u*(1.0-fburn)   + fburn*prob_parm.u_b;
  //vxt = 0.0; // initially at rest
  Real sumrhoY = 0.0;
  for (int n = 0; n < NUM_SPECIES; ++n) {
    Yt[n]   = Yu[n]*(1.0-fburn) + fburn*Yb[n];    
    sumrhoY += Yt[n];
  }
  // ensure sumY =1
  for (int n = 0; n < NUM_SPECIES; ++n) { Yt[n]   /=  sumrhoY;}
  // set T and P
  Tt = prob_parm.T_u*(1.0-fburn)   + fburn*prob_parm.T_b;
  Pt = prob_parm.p_u*(1.0-fburn)   + fburn*prob_parm.p_b ;
  
  // ignition spark
  Real Rspark = 5.0*prob_parm.lf;
  Real x0 =0; Real y0= prob_parm.Yflame;
  Real rad = sqrt((x- x0)*(x-x0) + (y- y0)*(y-y0));
  H = (rad-Rspark)/prob_parm.lf;
  const Real smoothspark=2.0;
  Real fspark = 0.5*tanh(smoothspark*H) + 0.5; // fspark   0:in  1:out
  Tt=1500.0*(1.0-fspark) + fspark*prob_parm.T_u;
  


  // compute density  
  cls.PYT2R(Pt,Yt,Tt,rhot);
  
  // compute energy  
  cls.RYP2E(rhot, Yt, Pt, et);

  //
  state(i, j, k, cls.UMX) = Real(0.0);
  state(i, j, k, cls.UMY) = rhot* vxt;
  state(i, j, k, cls.UMZ) = Real(0.0);
  state(i, j, k, cls.UET) = rhot * et + Real(0.5) * rhot * vxt * vxt;
  for (int n = 0; n < NUM_SPECIES; ++n) {
    state(i, j, k, cls.UFS + n) = rhot * Yt[n];
  }
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
  const int UFS  = ProbClosures::UFS;

  
  Real rhot, vxt,et,Tt,Pt;
  const Real *Yt;
   
  const int face = (idir+1)*sgn;

  switch(face)
  {
    case  2:  // SOUTH
      // inflow  unburn-----------------  
      rhot = prob_parm.rho_u;
      vxt  = prob_parm.u_u;
      Yt   = prob_parm.Y_u.data();
      et   = prob_parm.e_u;

      
     // Allow Pressure fluctuate  dP/dy = 0   T and Y FIX
      // rhot  = s_int[URHO]; vxt = s_int[UMY]/rhot;
      // et    = s_int[UET]/rhot - Real(0.5) * vxt* vxt/rhot;
      // closures.RYE2TP(rhot,Yt,et,Tt,Pt);   // give as P,T inside
      // Tt   = prob_parm.T_u;                // T fix
      // closures.PYT2R(Pt,Yt,Tt,rhot);       // recalculate rho
      // closures.RYP2E(rhot, Yt, Pt, et);    // recalvulate et

      // vxt = prob_parm.mflow/rhot;       // new vel


      s_ext[URHO] = rhot;
      s_ext[UMX]  = 0.0;
      s_ext[UMY]  = rhot* vxt;
      s_ext[UMZ]  = 0.0;    
      s_ext[UET]  = rhot * et + Real(0.5) * rhot * vxt * vxt;  
      for (int n = 0; n < NUM_SPECIES; ++n) {
        s_ext[UFS+n] = rhot * Yt[n];
      }        
      break;
    case  1:  // WEST      
      break;
    case -1:  // EAST
      break;
    case -2:   // NORTH
      break;
    default:

      break; 
  }
}
////////////////////////////////////////////////////////////////////////////////

// source term
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void user_source(
    int i, int j, int k, const auto &state, const auto &rhs,
    const ProbParm &lprobparm, ProbClosures const &closures, auto const dx) {}
////////////////////////////////////////////////////////////////////////////////

///////////////////////////////AMR//////////////////////////////////////////////
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void user_tagging(
    int i, int j, int k, int nt_level, auto &tagfab, const auto &sdatafab,
    const auto &geomdata, const ProbParm &prob_parm, int level) {
      
    const Real *prob_lo = geomdata.ProbLo();
    const Real *prob_hi = geomdata.ProbHi();
    const Real *dx = geomdata.CellSize();
    Real x = prob_lo[0] + (i + Real(0.5)) * dx[0];
    Real y = prob_lo[1] + (j + Real(0.5)) * dx[1];
   

    // SPECIES  (Li and Dryer mechanism)
    //              H2 O2 H2O H O OH HO2 H2O2 N2
    // index  UFS+  0   1  2  3 4  5  6   7   8

    
    Real Q[ProbClosures::NPRIM],U[ProbClosures::NCONS];
    for (int n = 0; n < ProbClosures::NCONS; ++n) {
      U[n] = sdatafab(i,j,k,n);
    }  

    auto thermo = ProbClosures::multispecies_pele_gas_t();    
    thermo.cons2prims_point(U,Q);

    Real rho = 0.0;
    for (int n = 0; n < NUM_SPECIES; ++n) {
          rho += sdatafab(i,j,k,ProbClosures::UFS + n);
    }  
    Real c =  (Q[ProbClosures::QT]-prob_parm.T_u)/(prob_parm.T_b - prob_parm.T_u);

    constexpr Real Cmax  = 0.95;constexpr Real Cmin  = 0.05;
    bool refine_flame =   (c < Cmax) && (c > Cmin);

    switch (level)
    {
      case 0:        
        tagfab(i,j,k) = refine_flame;
        break;
      case 1:
        tagfab(i,j,k) = refine_flame;
        break;
      default: 
        tagfab(i,j,k) = refine_flame;
        break;
    }

}
////////////////////////////////////////////////////////////////////////////////

}  // namespace PROB
#endif