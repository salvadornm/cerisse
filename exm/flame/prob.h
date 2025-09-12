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

#include <numbers>
#include <bc_types.h>


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
  Real rho_u = 0.9893277276003428;   // density  [kg/m^3]
  Real T_u   = 298;                   // temperature [K] 
  Real p_u   = 1.0132e+05;            // pressure [Pa]  (5 atm)
  GpuArray<Real, NUM_SPECIES> Y_u = { 0.014468, 0.22963 , 0.,0., 0., 0., 0.,0., 0.7559};  // mass fractions [-] 

  Real Y_0u[NUM_SPECIES] = {0.0};
  ProbParm(){
    // unburn  
    Y_0u[H2_ID]  = 0.014467517837899192;
    Y_0u[O2_ID]  = 0.22962878760502656;
    Y_0u[N2_ID]  = 0.7559036945570743;
  }



  // burn gases
  Real rho_b = 0.19626;                 // density  [kg/m^3]
  Real T_b   = 1644.8;                  // temperature [K]  
  Real p_b   = p_u;                     // pressure [Pa]  (1 atm)   
  GpuArray<Real, NUM_SPECIES> Y_b = {5.1557e-07,0.11471,0.12917,8.5022e-09,4.5034e-06,0.00021046,4.7965e-07,3.3189e-08,0.7559};

  // geometrical parameters                                     
  Real Lx     =   0.04;  // half-width domain
  Real Ly     =   0.04;
  Real Yflame =   0.5*Ly;

  Real pertur = 0.04; // perturbation of flame front based on flame thickness (0.1*lf)  

  Real SL     =  0.49; // estimated burning velocity
  Real lf     =  417e-6; // estimated flame thickness  (417 microns) Using Cantera and 1d-flame-plot.py
  // unburn gases velocity
  Real u_u     = 0.534292484155303; // inflow velocity (unburn) FROM CANTERA 1D profiles
    
  Real Q =  rho_u*u_u;  // incoming flow rate (per area)


};

// numerical method parameters
struct methodparm_t {

  public:

  static constexpr int  order = 2;        
  static constexpr Real use_LES = false;

};



using ProbClosures = closures_dt< indicies_t, transport_Pele_t, multispecies_pele_gas_t<indicies_t> >;

//using ProbRHS = rhs_dt< riemann_t<false, ProbClosures>, viscous_t<methodparm_t, ProbClosures>, reactor_t<ProbClosures> >;
//using ProbRHS = rhs_dt< riemann_t<false, ProbClosures>, viscous_t<methodparm_t, ProbClosures>, no_source_t>;
//using ProbRHS = rhs_dt< riemann_t<false, ProbClosures>, no_diffusive_t, reactor_t<ProbClosures> >;
//using ProbRHS = rhs_dt< no_euler_t, viscous_t<methodparm_t, ProbClosures>, reactor_t<ProbClosures> >;
using ProbRHS = rhs_dt< riemann_t<false, ProbClosures>, viscous_t<methodparm_t, ProbClosures>, reactor_t<ProbClosures> >;

//using ProbRHS = rhs_dt< centraldif_t <false, false, 2,ProbClosures>, viscous_t<methodparm_t, ProbClosures>, reactor_t<ProbClosures> >;



typedef manual_bc_t<ProbClosures> GlobalBC;


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
    ProbParm const &prob_parm, Utility* util = nullptr) {
  const Real *prob_lo = geomdata.ProbLo();
  const Real *prob_hi = geomdata.ProbHi();
  const Real *dx = geomdata.CellSize();

  Real x = prob_lo[0] + (i + Real(0.5)) * dx[0];
  Real y = prob_lo[1] + (j + Real(0.5)) * dx[1];
  
  Real rhot, vxt,et,Tt,Pt;
  Real Yt[NUM_SPECIES];const Real *Yu;const Real *Yb;


  Real yinterf =  prob_parm.Yflame;
  
  // smooth profile over lf 
  // const Real smooth=5; // the smaller the wider the spread (smooth 5, spread over -lf:lf)
  // Real H = (y-yinterf)/prob_parm.lf;
  // Real fburn = 0.5*tanh(smooth*H) + 0.5; // fburn   0:unburn  1:burn

  // Constant Pressure 
  Pt = prob_parm.p_u ;
  
  //--- Initialise 1D flame from data -------------------------------------------------------------------
  // PMF data from 0...to 0.04 (Flame front position: 0.015634 m) obtained from tools/combustion/1d-flame-plot.py
  Real yflame_front = 0.015634;
  // upper and lower position of the cell (relative to flame)
  Real y1  = y - 0.5*dx[1]- yinterf;
  Real y2  = y + 0.5*dx[1]- yinterf;
  // shift to flame_front
  y1 += yflame_front;
  y2 += yflame_front;

  Real sumrhoY = 0.0;
  
  // read from PMF profile ----> pmf_vals
  //--------------------------------------------------------------------------------------
  // pmf_vals[0] =T  pmf_vals[1]= Velocity  pmf_vals[2] = rho pmf_vals[3+k] = Y[k];
  GpuArray<Real, NUM_SPECIES + 4 > pmf_vals = {0.0}; 

  pele::physics::PMF::PmfData::DataContainer *pmf_data = util->pmfData.getDeviceData();
  pele::physics::PMF::pmf(pmf_data,y1,y2,pmf_vals);

  // PMF--> T,u and Y (P is assumed constant)
  // Tt  = pmf_vals[0];
  // vxt = pmf_vals[1];
  // for (int n = 0; n < NUM_SPECIES; ++n) {
  //   Yt[n]   = max(pmf_vals[3+n],0.0);        
  //   sumrhoY += Yt[n];
  // }
  
  
  // fresh start
  Tt  = prob_parm.T_u;
  vxt = prob_parm.u_u;
  for (int n = 0; n < NUM_SPECIES; ++n) {
    Yt[n] = prob_parm.Y_0u[n];
    sumrhoY += Yt[n];
  }
  // spark
  // Gaussian profile centered at interface, width 0.005
  const amrex::Real center = yinterf;
  const amrex::Real sigma  = 0.005;
  // Gaussian value between 0 and 1
  amrex::Real gauss = std::exp( - std::pow((y - center)/sigma, 2) );
  Tt += 1000.0*gauss;



  // ensure sumY =1
  for (int n = 0; n < NUM_SPECIES; ++n) { Yt[n]   /=  sumrhoY;}
  //--------------------------------------------------------------------------------------
  
  // compute density  
  cls.PYT2R(Pt,Yt,Tt,rhot);
  
  // compute energy  
  cls.RYP2E(rhot, Yt, Pt, et);

  // debug
  Real haux[NUM_SPECIES]={0.0};
  cls.RTY2Hi(rhot, Tt, Yt, haux);
  Real ht = 0.0; 
  Real kin = Real(0.5) * vxt * vxt;

  // std::cout << " ------------------------------- " << std::endl;
  // std::cout << " (spec internal energy) e "  << et << " [J/kg] " << std::endl;
  // std::cout << " T "    << Tt << " [kg/m3] " << std::endl;
  // std::cout << " rho "  << rhot << " [kg/m3] " << std::endl;
  // std::cout << " 1/2 v^2     "  << kin << " [J/kg] " << std::endl;
  // std::cout << " 1/2 rho v^2 "  << rhot*kin << " [J/m3] " << std::endl;
  // std::cout << " spec total Energy=  et " << et + kin << " [J/m3] " << std::endl;
  // std::cout << " total Energy= rho et " << rhot * (et + kin )<< " [J/m3] " << std::endl;
  // for (int n = 0; n < NUM_SPECIES; ++n) {
  //  std::cout << " n "  <<  n << " hk = " << haux[n] << " [J/kg] " << std::endl;
  //  ht += haux[n]*Yt[n];
  // }
  // std::cout << " (spec enthalpy M1) h  = Yk hk      => "  << ht << " [J/kg] " << std::endl;
  // std::cout << " (spec enthalpy M2) h  = e + P/rho  => "  << et + Pt/rhot<< " [J/kg] " << std::endl;
  // std::cout << " (spec total enthalpy) htot  = et + P/rho  => "  << et + Pt/rhot + kin << " [J/kg] " << std::endl;

  // // fix total enthalpy
  // Real htotalfix = -130.2384975 + Real(0.5) * prob_parm.u_u* prob_parm.u_u;
  // Real hfix= htotalfix - kin;
  // Real Tfix = 0.0;
  // // calc new T  and rho
  // auto eos = pele::physics::PhysicsType::eos();
  // eos.RHY2T(rhot*rho_si2cgs,hfix*specenergy_si2cgs, Yt, Tfix); cls.PYT2R(Pt,Yt,Tfix,rhot);
  // // re-calculate internal energy
  // Real efix = hfix - Pt/rhot;
  // et = efix;  
  /// SNM

  //
  state(i, j, k, cls.UMX) = Real(0.0);
  state(i, j, k, cls.UMY) = rhot* vxt;
  state(i, j, k, cls.UMZ) = Real(0.0);
  state(i, j, k, cls.UET) = rhot * et + Real(0.5) * rhot * vxt * vxt;
  for (int n = 0; n < NUM_SPECIES; ++n) {
    state(i, j, k, cls.UFS + n) = rhot * Yt[n];
  }


  //

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

      //GlobalBC::bc_fixP(0.0,1.0,0.0,&closures,prob_parm.p_u, s_int, s_ext);

      break;
    case  1:  // WEST      
      break;
    case -1:  // EAST
      break;      
    case -2:   // NORTH
    //  GlobalBC::bc_fixP(0.0,-1.0,0.0,&closures,prob_parm.p_u, s_int, s_ext);
      GlobalBC::bc_subsonic_outflow_fixP(0.0,-1.0,0.0,&closures,prob_parm.p_u, s_int, s_ext);
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
