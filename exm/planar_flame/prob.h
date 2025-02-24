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
//       temperature   700 K
//          pressure   5.0662e+05 Pa
//           density   2.1059 kg/m^3
//  mean mol. weight   24.192 kg/kmol
//
//                         1 kg             1 kmol     
//                     ---------------   ---------------
//          enthalpy        4.9633e+05        1.2007e+07  J
//   internal energy        2.5575e+05        6.1872e+06  J
//                      mass frac. Y      mole frac. X     chem. pot. / RT
//                     ---------------   ---------------   ---------------
//                H2          0.014468           0.17361           -16.827
//                O2           0.22963           0.17361           -25.822
//                N2            0.7559           0.65278           -22.838


// Burning composition obtained with Cantera equivalence ratio 0.5
// using  cerisse/tools/combustion/equil_from_equivalenceratio.py
// 
//       temperature   1980.5 K
//          pressure   5.0663e+05 Pa
//           density   0.81467 kg/m^3
//  mean mol. weight   26.479 kg/kmol
//
//                          1 kg             1 kmol     
//                     ---------------   ---------------
//          enthalpy        4.9633e+05        1.3142e+07  J
//   internal energy       -1.2555e+05       -3.3243e+06  J
//                      mass frac. Y      mole frac. X     chem. pot. / RT
//                    ---------------   ---------------   ---------------
//                H2         5.172e-06        6.7932e-05           -27.423
//                O2           0.11432          0.094605           -29.468
//               H2O           0.12869           0.18915           -42.156
//                 H        1.9846e-07        5.2134e-06           -13.711
//                 O        4.7351e-05        7.8369e-05           -14.734
//                OH         0.0010298         0.0016034           -28.445
//               HO2        3.7555e-06        3.0129e-06           -43.179
//              H2O2        2.9303e-07        2.2812e-07            -56.89
//                N2            0.7559           0.71449           -25.625

// problem parameters
struct ProbParm {

  // unburn gases
  Real rho_u = 2.1059;               // density  [kg/m^3]
  Real T_u   =  700;                 // temperature [K] 
  Real u_u   = 5.3;                  // velocity [m/s]
  Real p_u   = 5.0663e+05;           // pressure [Pa]  (5 atm)
  Real e_u   = 2.5575e+05;         // internal energy [J/kg]  
 
  combustion_functions::Y_Mixture massfr = combustion_functions::mol2Y_H2Air(mix);  // unburn mass fractions
  GpuArray<Real, NUM_SPECIES> Y_u = { massfr.H2, massfr.O2, 0.,0., 0., 0., 0.,0., massfr.N2};  // mass fractions [-] 

  // burn gases
  Real rho_b = 0.81467;                 // density  [kg/m^3]
  Real T_b = 1980.5;                    // temperature [K]  
  Real u_b = 0.0;                       // velocity [m/s]
  Real p_b = p_u;                       // pressure [Pa]  (1 atm)   
  Real e_b = -1.2555e+05 ;              // internal energy [J/kg]
 
  // butn mass fractions
  GpuArray<Real, NUM_SPECIES> Y_b = {5.172e-06,0.11432,0.12869,1.9846e-07,4.7351e-05,0.0010298,3.7555e-06,2.9303e-07,0.7559};

  // geometrical parameters                                     
  Real Yflame = 0.05;
  Real lf = 507e-6; // flame thicklness  (507 microns)
  Real Ly = 0.025;  // width domain


};

using ProbClosures = closures_dt< indicies_t, visc_suth_t, cond_suth_t, multispecies_perfect_gas_t<indicies_t> >;
//using ProbRHS = rhs_dt< weno_t<ReconScheme::Teno5, ProbClosures>, no_diffusive_t, reactor_t<ProbClosures> >;
using ProbRHS = rhs_dt< riemann_t<false, ProbClosures>, no_diffusive_t, reactor_t<ProbClosures> >;


void inline inputs() {
  // ParmParse pp;  
  amrex::Print() << " ****** Starting  *******" <<  std::endl;
  amrex::Print() << " Planar H2-Air Flame (February 2025)  " <<  std::endl;
  amrex::Print() << " Nspecies = " << NUM_SPECIES << std::endl;
  amrex::Print() << " H2 O2 H2O H O OH HO2 H2O2 N2 " << std::endl;
  amrex::Print() << " ******           *******" << std::endl;

  ProbParm const prob_parm;
  const Real *Yt;
  
  amrex::Print() << " composition unburn " << std::endl;
  amrex::Print() << " ------------------ " << std::endl;

  Yt = prob_parm.Y_u.data();
  amrex::Real sumY =0.0;
  for (int n = 0; n < NUM_SPECIES; ++n) {    
    amrex::Print() << "n= " <<  n << " Y= "<< Yt[n]<< std::endl;
    sumY += Yt[n];
  }
  amrex::Print() << " ------------------ " << std::endl;
  amrex::Print() << " sum Y=  " << sumY << std::endl;
  
  amrex::Print() << " composition burn " << std::endl;
  amrex::Print() << " ------------------ " << std::endl;
  
  Yt = prob_parm.Y_b.data();
  amrex::Real sumY2 =0.0;
  for (int n = 0; n < NUM_SPECIES; ++n) {    
    amrex::Print() << "n= " <<  n << " Y= "<< Yt[n]<< std::endl;
    sumY2 += Yt[n];
  }
  amrex::Print() << " ------------------ " << std::endl;
  amrex::Print() << " sum Y=  " << sumY2 << std::endl;

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
  
  Real rhot, vxt,et;
  const Real *Yt;

  Real wavelength = prob_parm.Ly;
  Real wave= 2.0*std::numbers::pi/wavelength;
  // pertubastion and frequency
  const Real ypertur = 0.005;//4*dx[1]; 
  const Real Nwaves = 1;

  Real yinterf =  prob_parm.Yflame  + ypertur*cos(Nwaves*wave*x);
  
  // compute 

  if (y < yinterf) {   // unburn mixture at the bottom    
    rhot = prob_parm.rho_u;
    vxt  = prob_parm.u_u;
    Yt   = prob_parm.Y_u.data();
    et   = prob_parm.e_u;
  } else {            // burn mixture at the top    
    rhot = prob_parm.rho_b;
    vxt  = prob_parm.u_b;
    Yt   = prob_parm.Y_b.data();
    et   = prob_parm.e_b;
  }
  
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
  
  Real rhot, vxt,et;
  const Real *Yt;
   
  const int face = (idir+1)*sgn;

  switch(face)
  {
    case  2:  // SOUTH
      // inflow  unburn      
      rhot = prob_parm.rho_u;
      vxt  =  prob_parm.u_u;
      Yt   = prob_parm.Y_u.data();
      et   = prob_parm.e_u;
      //Pt   = prob_parm.p_u;
      //closures.RYP2E(rhot, Yt, Pt, et);

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
  
    constexpr Real threshold = 0.1;

    // tag cells based on gradients of H2 (change here)
    const int QG = ProbClosures::UFS+0;
    
    Real T = sdatafab(i,j,k,QG) + max(1e-5*sdatafab(i,j,k,QG),1e-7);
    Real dTx = Math::abs( sdatafab(i+1,j,k,QG)- sdatafab(i,j,k,QG)   )/T;
    Real dTy = Math::abs( sdatafab(i,j+1,k,QG)- sdatafab(i,j-1,k,QG) )/T;
    Real gradT = sqrt(dTx*dTx + dTy*dTy);

    switch (level)
    {
      case 0:
        tagfab(i,j,k) = (gradT > threshold);
        break;
      case 1:
        tagfab(i,j,k) = (gradT > threshold);
        break;
      default:
        tagfab(i,j,k) = (gradT > threshold);
        break;
    }

}
////////////////////////////////////////////////////////////////////////////////

}  // namespace PROB
#endif