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


// SPECIES  (Li and Dryer order)
// H2 O2 H2O H O OH HO2 H2O2 N2

// problem parameters
struct ProbParm {

  // unburn gases
  Real rho_u = 1.2;                  // density  [kg/m^3]
  Real T_u   =  700;                 // temperature [K] 
  Real u_u   = 5.3;                  // velocity [m/s]
  Real p_u   = 1.0*pres_atm2si;      // pressure [Pa]  (1 atm)

  combustion_functions::Y_Mixture massfr = combustion_functions::mol2Y_H2Air(mix);  // unburn mass fractions
  GpuArray<Real, NUM_SPECIES> Y_u = { massfr.H2, massfr.O2, 0.,0., 0., 0., 0.,0., massfr.N2};  // mass fractions [-] 

  // burn gases
  Real rho_b = 0.18075;               // density  [kg/m^3]
  Real T_b = 1973;                    // temperature [K]  
  Real u_b = 0.0;                     // velocity [m/s]
  Real p_b = p_u;                     // pressure [Pa]  (1 atm)   

  GpuArray<Real, NUM_SPECIES> Y_b = {0.,         0.01277243, 0., 0., 0., 0.10136214, 0.,         0., 0.};

  // geometrical parameters                                     
  Real Yflame = 0.05;
  Real lf = 507e-6; // flame thicklness  (507 microns)


};

using ProbClosures = closures_dt< indicies_t, visc_suth_t, cond_suth_t, multispecies_perfect_gas_t<indicies_t> >;
using ProbRHS = rhs_dt< weno_t<ReconScheme::Teno5, ProbClosures>, no_diffusive_t, reactor_t<ProbClosures> >;

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
  
  Real Pt, rhot, vxt;
  const Real *Yt;

  Real wavelength = 10*prob_parm.lf;
  Real wave= 2.0*std::numbers::pi/wavelength;
  Real ypertur = 2*dx[1];
  Real yinterf =  prob_parm.Yflame;  + ypertur*cos(wave*x);


  // compute 


  if (y < yinterf) {   // unburn data
    Pt = prob_parm.p_u;
    rhot = prob_parm.rho_u;
    vxt = prob_parm.u_u;
    Yt = prob_parm.Y_u.data();
  } else {
    Pt = prob_parm.p_b;
    rhot = prob_parm.rho_b;
    vxt = prob_parm.u_b;
    Yt = prob_parm.Y_b.data();
  }
  Real et;
  cls.RYP2E(rhot, Yt, Pt, et);

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
  
  Real Pt, rhot, vxt,et;
  const Real *Yt;
   
  const int face = (idir+1)*sgn;

  switch(face)
  {
    case  2:  // SOUTH
      // inflow
      Pt   = prob_parm.p_u;
      rhot = prob_parm.rho_u;
      vxt  =  prob_parm.u_u;
      Yt   = prob_parm.Y_u.data();
      closures.RYP2E(rhot, Yt, Pt, et);

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
  // Real dengrad_threshold = 0.5;
  // Real drhox = Math::abs(sdatafab(i+1,j,k,URHO) -
  // sdatafab(i-1,j,k,URHO))/sdatafab(i,j,k,URHO); if (drhox >
  // dengrad_threshold) {
  //   tagfab(i,j,k) = true;
  //   tagfab(i+1,j,k) = true;
  //   tagfab(i+2,j,k) = true;
  //   tagfab(i+3,j,k) = true;
  // }
}
////////////////////////////////////////////////////////////////////////////////

}  // namespace PROB
#endif