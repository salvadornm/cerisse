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

// executing python equil_from_equivalenceratio.py  

//  Initial Conditions  Mixture
//     P=  101325.0 [Pa] amd T =  298.0 [K]
//     Equivalence Ratio=  0.8
//  Unburn Mixture ...
// ************ Phase gas ************
// Moles:  1.0

//   gas:

//        temperature   298 K
//           pressure   1.0132e+05 Pa
//            density   0.90377 kg/m^3
//   mean mol. weight   22.1 kg/kmol
//    phase of matter   gas

//                           1 kg             1 kmol     
//                      ---------------   ---------------
//           enthalpy           -137.15             -3031  J
//    internal energy       -1.1225e+05       -2.4807e+06  J
//            entropy            8424.8        1.8619e+05  J/K
//     Gibbs function       -2.5107e+06       -5.5487e+07  J
//  heat capacity c_p            1314.9             29059  J/K
//  heat capacity c_v            938.65             20744  J/K

//                       mass frac. Y      mole frac. X     chem. pot. / RT
//                      ---------------   ---------------   ---------------
//                 H2          0.022949           0.25157           -17.086
//                 O2           0.22765           0.15723           -26.512
//                 N2            0.7494           0.59119           -23.559
//      [   +6 minor]                 0                 0  

//  Burn Mixture ...
// ************ Phase gas ************
// Moles:  0.8763651682839297

//   gas:

//        temperature   2175.1 K
//           pressure   1.0133e+05 Pa
//            density   0.14129 kg/m^3
//   mean mol. weight   25.218 kg/kmol
//    phase of matter   gas

//                           1 kg             1 kmol     
//                      ---------------   ---------------
//           enthalpy           -137.15           -3458.6  J
//    internal energy       -7.1729e+05       -1.8089e+07  J
//            entropy             10544         2.659e+05  J/K
//     Gibbs function       -2.2935e+07       -5.7836e+08  J
//  heat capacity c_p            1619.1             40830  J/K
//  heat capacity c_v            1289.4             32515  J/K

//                       mass frac. Y      mole frac. X     chem. pot. / RT
//                      ---------------   ---------------   ---------------
//                 H2         0.0001171         0.0014648           -26.263
//                 O2          0.044507          0.035076           -32.406
//                H2O             0.202           0.28276           -42.466
//                  H         7.435e-06        0.00018601           -13.131
//                  O        0.00027133        0.00042768           -16.203
//                 OH         0.0036952         0.0054792           -29.335
//                HO2        3.1809e-06        2.4303e-06           -45.538
//               H2O2        2.3096e-07        1.7123e-07           -58.669
//                 N2            0.7494            0.6746           -27.612
// problem parameters
struct ProbParm {

  // unburn gases
  Real rho_u = 0.90377;               // density  [kg/m^3]
  Real T_u   = 298;                   // temperature [K] 
  Real p_u   = 1.0132e+05;            // pressure [Pa]  (5 atm)
  Real e_u   = -1.1225e+05;           // internal energy [J/kg]  

  // burn gases
  Real rho_b = 0.14129;                 // density  [kg/m^3]
  Real T_b   = 2175.1;                  // temperature [K]  
  Real p_b   = p_u;                     // pressure [Pa]  (1 atm)   
  Real e_b   = -7.1729e+05 ;            // internal energy [J/kg]
  GpuArray<Real, NUM_SPECIES> Y_b = {0.0001171,0.044507,0.202 ,7.435e-06,0.00027133,0.0036952,3.1809e-06,2.3096e-07,0.7559};

  Real Y_0[NUM_SPECIES] = {0.0};
  ProbParm(){
    Y_0[H2_ID] = 0.022949 ;
    Y_0[O2_ID] = 0.22765;
    Y_0[N2_ID] = 1.0 - Y_0[O2_ID] - Y_0[H2_ID];
  }
  
  GpuArray<Real, NUM_SPECIES> Y_u = { Y_0[H2_ID] ,Y_0[O2_ID] , 0.,0., 0., 0., 0.,0., Y_0[N2_ID]};  // mass fractions [-] 


  // geometrical parameters                                     
  Real Lx     =   0.04;  // half-width domain
  Real Ly     =   0.04;
  Real Yflame =   0.5*Ly;

  // unburn gases velocity
  Real u_u     = 3.2; // inflow velocity (unburn) IMPOSED [m/s]

  
  Real Q =  rho_u*u_u;  // flow rate (per area)

  
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

using ProbRHS = rhs_dt< riemann_t<false, ProbClosures>, viscous_t<methodparm_t, ProbClosures>, reactor_t<ProbClosures> >;
//using ProbRHS = rhs_dt< riemann_t<false, ProbClosures>, viscous_t<methodparm_t, ProbClosures>, no_source_t >;
//using ProbRHS = rhs_dt< no_euler_t, viscous_t<methodparm_t, ProbClosures>, no_source_t >;


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
  amrex::Print() << " Planar Bunsen H2-Air Flame (June 2025)  " <<  std::endl;
  amrex::Print() << " Nspecies = " << NUM_SPECIES << std::endl;
  amrex::Print() << " H2 O2 H2O H O OH HO2 H2O2 N2 " << std::endl;
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
  Real yinterf= -0.1/1000.0; // [m]

  
  Yu  = prob_parm.Y_u.data();
  Yb  = prob_parm.Y_b.data();

  vxt = prob_parm.u_u;
  Real sumrhoY = 0.0;

  if (y > yinterf)
  { // burn
    Tt = prob_parm.T_b;
    for (int n = 0; n < NUM_SPECIES; ++n) {
      Yt[n]    = Yb[n];    
      sumrhoY += Yt[n];
    }
  }
  else
  { // unburn
    Tt  = prob_parm.T_u;
    for (int n = 0; n < NUM_SPECIES; ++n) {
      Yt[n]   = prob_parm.Y_0[n];    
      sumrhoY += Yt[n];
    }

  }

  

  // ensure sumY =1
  for (int n = 0; n < NUM_SPECIES; ++n) { Yt[n]   /=  sumrhoY;}
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
        prob_parm.Q,prob_parm.T_u,prob_parm.Y_0, s_int, s_ext);            
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
