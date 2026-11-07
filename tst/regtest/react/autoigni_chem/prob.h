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
static combustion_functions::Moles_Mixture mix = combustion_functions::compute_MixMol_fromEquivalenceRatio_H2Air(Phi);

// SPECIES  (Li and Dryer mechanism)
// H2 O2 H2O H O OH HO2 H2O2 N2
// Fresh gases
//
//                     mass frac. Y      mole frac. X     
//                     ---------------   ---------------   
//                H2          0.014468           0.17361  
//                O2           0.22963           0.17361  
//                N2            0.7559           0.65278


// problem parameters
struct ProbParm {
  // unburn gases
  Real T0   = 1000;                 // temperature [K] 
  Real p0   = 1.0132e+05;           // pressure [Pa]  (1 atm)
  combustion_functions::Y_Mixture massfr = combustion_functions::mol2Y_H2Air(mix);  // unburn mass fractions
  GpuArray<Real, NUM_SPECIES> Y0 = { 0.014468, 0.22963 , 0.,0., 0., 0., 0.,0., 0.7559};  // mass fractions [-] 
};

using ProbClosures = closures_dt< indicies_t, transport_Pele_t, multispecies_pele_gas_t<indicies_t> >;
using ProbRHS = rhs_dt< no_euler_t, no_diffusive_t, reactor_t<ProbClosures> >;

void inline inputs() {
  // ParmParse pp;  
  amrex::Print() << " ****** Starting  *******" <<  std::endl;
  amrex::Print() << "  H2-Air Flame (Autoignition)  " <<  std::endl;
  amrex::Print() << " ******           *******" << std::endl;
}

// initial condition
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void prob_initdata(
    int i, int j, int k, Array4<Real> const &state,
    GeometryData const &geomdata, ProbClosures const &cls,
    ProbParm const &prob_parm) {
    
  Real rhot, et,Tt,Pt;
  const Real *Yt;

  // initial conditions
  Yt = prob_parm.Y0.data();
  Tt = prob_parm.T0;
  Pt = prob_parm.p0; 

  // compute density  
  cls.PYT2R(Pt,Yt,Tt,rhot);  
  // compute energy  
  cls.RYP2E(rhot, Yt, Pt, et);

  state(i, j, k, cls.UMX) = Real(0.0);
  state(i, j, k, cls.UMY) = Real(0.0);
  state(i, j, k, cls.UMZ) = Real(0.0);
  state(i, j, k, cls.UET) = rhot * et;
  for (int n = 0; n < NUM_SPECIES; ++n) {
    state(i, j, k, cls.UFS + n) = rhot * Yt[n];
  }
}

/////////////////////////////// BC /////////////////////////////////////////////
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
bcnormal(const Real x[AMREX_SPACEDIM], Real dratio, const Real s_int[ProbClosures::NCONS],
         const Real s_refl[ProbClosures::NCONS], Real s_ext[ProbClosures::NCONS], const int idir,
         const int sgn, const Real time, GeometryData const & /*geomdata*/,
         ProbClosures const &closures, ProbParm const &prob_parm) {}
////////////////////////////////////////////////////////////////////////////////

// source term
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void user_source(
    int i, int j, int k, const auto &state, const auto &rhs,
    const ProbParm &lprobparm, ProbClosures const &closures, auto const dx) {}
////////////////////////////////////////////////////////////////////////////////

///////////////////////////////AMR//////////////////////////////////////////////
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void user_tagging(
    int i, int j, int k, int nt_level, auto &tagfab, const auto &sdatafab,
    const auto &geomdata, const ProbParm &prob_parm, int level) {}
////////////////////////////////////////////////////////////////////////////////

}  // namespace PROB
#endif