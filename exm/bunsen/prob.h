#ifndef CNS_PROB_H_
#define CNS_PROB_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_ParmParse.H>
#include <PelePhysics.H>
#include <ReactorBase.H>

#include "Constants.h"
#include "CombustionFunctions.h"
#include "Utilities.h"

#include "Closures.h"
#include "RHS.h"

#if CNS_USE_EB     
#include <ebm.h>
#include <walltypes.h>
#endif
  
#include <bc_types.h>



namespace PROB {

// executing python equil_from_equivalenceratio.py (ubnder Kee/) 

//  Initial Conditions  Mixture
//     P=  101325.0 [Pa] amd T =  300.0 [K]
//     Equivalence Ratio=  0.8
///    Fuel: 87.5% H2  12.5% CH4 (per volume)


//  Unburn Mixture ...
// ************ Phase gas ************
// Moles:  1.0

//   gas:

//        temperature   300 K
//           pressure   1.0133e+05 Pa
//            density   0.97189 kg/m^3
//   mean mol. weight   23.925 kg/kmol
//    phase of matter   gas

//                           1 kg             1 kmol     
//                      ---------------   ---------------
//           enthalpy            -74239       -1.7762e+06  J
//    internal energy       -1.7849e+05       -4.2705e+06  J
//            entropy            8007.6        1.9158e+05  J/K
//     Gibbs function       -2.4765e+06       -5.9251e+07  J
//  heat capacity c_p            1222.7             29254  J/K
//  heat capacity c_v            875.19             20939  J/K

//                       mass frac. Y      mole frac. X     chem. pot. / RT
//                      ---------------   ---------------   ---------------
//                CH4           0.01646          0.024547            -56.03
//                 H2          0.014479           0.17183           -17.479
//                 O2            0.2257           0.16876           -26.453
//                 N2           0.74336           0.63486           -23.487

//  Burn Mixture ...
// ************ Phase gas ************
// Moles:  0.915684304383077

//   gas:

//        temperature   2117.2 K
//           pressure   1.0133e+05 Pa
//            density   0.15039 kg/m^3
//   mean mol. weight   26.128 kg/kmol
//    phase of matter   gas

//                           1 kg             1 kmol     
//                      ---------------   ---------------
//           enthalpy            -74239       -1.9397e+06  J
//    internal energy       -7.4797e+05       -1.9543e+07  J
//            entropy             10197        2.6643e+05  J/K
//     Gibbs function       -2.1663e+07       -5.6602e+08  J
//  heat capacity c_p            1561.4             40796  J/K
//  heat capacity c_v            1243.2             32482  J/K

//                       mass frac. Y      mole frac. X     chem. pot. / RT
//                      ---------------   ---------------   ---------------
//               CH2O         1.796e-13        1.5628e-13           -68.619
//                HCO        1.4265e-11        1.2844e-11           -55.242
//                CO2           0.04438          0.026348           -58.005
//                 CO        0.00049205          0.000459           -41.865
//                 H2        6.4213e-05        0.00083223           -26.753
//                  H        3.8577e-06        9.9994e-05           -13.377
//                 O2          0.044704          0.036503            -32.28
//                  O        0.00018305        0.00029895            -16.14
//                 OH         0.0022825         0.0035066           -29.517
//                HO2        2.3351e-06        1.8485e-06           -45.657
//               H2O2        1.6022e-07        1.2307e-07           -59.033
//                H2O           0.16453           0.23863           -42.893
//                 N2           0.74336           0.69332           -27.491
//      [   +4 minor]        5.9832e-20        1.0342e-19  
//--------------------------------------------------------------------------
// problem parameters
struct ProbParm {

  // unburn gases
  static constexpr Real rho_u = 0.97189;               // density  [kg/m^3]
  static constexpr Real T_u  = 300.0;                   // temperature [K] 
  static constexpr Real p_u   = 1.01325e+05;            // pressure [Pa]  (5 atm)

  // burn gases
  static constexpr Real rho_b = 0.15039 ;                 // density  [kg/m^3]
  static constexpr Real T_b   = 2117.2;                  // temperature [K]  
  static constexpr Real p_b  = p_u;                     // pressure [Pa]  (1 atm)   

  // mass fractions
  
  Real Y_0b[NUM_SPECIES]   = {0.0};
  Real Y_0u[NUM_SPECIES] = {0.0};

  ProbParm(){
    // unburn
    Y_0u[CH4_ID] = 0.01646;
    Y_0u[H2_ID]  = 0.014479;
    Y_0u[O2_ID]  = 0.2257;
    Y_0u[N2_ID]  = 0.74336; 
    // burn
    Y_0b[CH2O_ID] = 1.796e-13;
    Y_0b[HCO_ID]  = 1.4265e-11;
    Y_0b[CO2_ID]  = 0.04438;
    Y_0b[CO_ID]   = 0.00049205;
    Y_0b[H2_ID]   = 6.4213e-05;
    Y_0b[H_ID]    = 3.8577e-06;
    Y_0b[O2_ID]   = 0.044704;
    Y_0b[O_ID]    = 0.00018305;
    Y_0b[OH_ID]   = 0.0022825;
    Y_0b[HO2_ID]  = 2.3351e-06;
    Y_0b[H2O2_ID] = 1.6022e-07;
    Y_0b[H2O_ID]  = 0.16453;
    Y_0b[N2_ID]   = 0.74336;
  }
  
  // Initial Flame position                                
  Real Yflame= 0.8/1000.0; // [m] <-------------

  // unburn gases velocity [m/s] (imposed from Fruzza et al 2023)
  // Real u_u    = 3.2; // Flame with 100% H2
  Real u_u     = 2.0; // Flame with 87.5% H2  12.5% CH4

  Real u_b    = 6.7852789; // from 1D solution

  Real Q =  rho_u*u_u;  // incoming flow rate (per area)
  
};

// numerical method parameters
struct methodparm_t {

  public:

  static constexpr int  order = 2;        
  static constexpr Real use_LES = false;

};

struct wall_param {

  public:

  static constexpr Real Twall = 800.0;           // wall temperature (if isothermal used)
  static constexpr bool solve_diffwall = true;   // solve viscous effects at walls
  
};


//template <typename cls_t > class user_source_t;


using ProbClosures = closures_dt< indicies_t, transport_Pele_t, multispecies_pele_gas_t<indicies_t> >;

//using ProbRHS = rhs_dt< riemann_t<false, ProbClosures>, viscous_t<methodparm_t, ProbClosures>, reactor_t<ProbClosures> >;
//using ProbRHS = rhs_dt< riemann_t<false, ProbClosures>, no_diffusive_t, no_source_t >;
using ProbRHS = rhs_dt< riemann_t<false, ProbClosures>, viscous_t<methodparm_t, ProbClosures>, no_source_t >;


// define type of wall and EBM class
#if CNS_USE_EB    

//typedef adiabatic_wall_t<ProbClosures> TypeWall;

typedef isothermal_wall_t<wall_param,ProbClosures> TypeWall;

typedef ebm_t<TypeWall,wall_param,ProbClosures> ProbEB;
#endif

typedef manual_bc_t<ProbClosures> GlobalBC;


void inline inputs() {
  // ParmParse pp;  
  amrex::Print() << " ****** Starting  *******" <<  std::endl;
  amrex::Print() << " Planar Bunsen CH4/H2-Air Flame (June 2025)  " <<  std::endl;
  amrex::Print() << " Nspecies = " << NUM_SPECIES << std::endl;
  amrex::Print() << " ******           *******" << std::endl;
}

// initial condition
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void prob_initdata(
    int i, int j, int k, Array4<Real> const &state,
    GeometryData const &geomdata, ProbClosures const &cls,
    ProbParm const &prob_parm, Utility* util = nullptr) {
  const Real *prob_lo = geomdata.ProbLo();
  const Real *prob_hi = geomdata.ProbHi();
  const Real *dx = geomdata.CellSize();

  Real x = prob_lo[0] + (i + Real(0.5)) * dx[0];
  Real y = prob_lo[1] + (j + Real(0.5)) * dx[1];
  
  Real rhot, vxt,et,Tt,Pt;
  Real Yt[NUM_SPECIES];const Real *Yu;const Real *Yb;

 
  // Constant Pressure 
  Pt = prob_parm.p_u ;
  
  // initial position flame (domain coordinates)
  Real yinterf= prob_parm.Yflame; // [m] <-------------

  Real sumY = 0.0;

  if (y > yinterf)
  { // burn
    vxt = prob_parm.u_b;
    Tt  = prob_parm.T_b;   
    for (int n = 0; n < NUM_SPECIES; ++n) {
      Yt[n]    = prob_parm.Y_0b[n];          
      sumY += Yt[n];
    }
  }
  else
  { // unburn
    vxt = prob_parm.u_u;
    Tt  = prob_parm.T_u;
    for (int n = 0; n < NUM_SPECIES; ++n) {
      Yt[n]   = prob_parm.Y_0u[n];    
      sumY += Yt[n];
    }

  }


  // ensure sumY =1
  for (int n = 0; n < NUM_SPECIES; ++n) { Yt[n]   /=  sumY;}
  //--------------------------------------------------------------------------------------
  
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

  const int face = (idir+1)*sgn;

  switch(face)
  {
    case  2:  // SOUTH
      // inflow  unburn-----------------        
      GlobalBC::bc_inlet_fixmassflow(0.0,1.0,0.0,&closures,
        prob_parm.Q,prob_parm.T_u,prob_parm.Y_0u, s_int, s_ext);            
      break;
    case  1:  // WEST      
      break;
    case -1:  // EAST
      break;
    case -2:   // NORTH (fix pressure)
      GlobalBC::bc_fixP(0.0,-1.0,0.0,&closures,prob_parm.p_u, s_int, s_ext);
      break;
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
    //Real x = prob_lo[0] + (i + Real(0.5)) * dx[0];
    //Real y = prob_lo[1] + (j + Real(0.5)) * dx[1];
   

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

    constexpr Real Cmax  = 0.9;constexpr Real Cmin  = 0.1;
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
