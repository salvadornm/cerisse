#ifndef PASR_H_
#define PASR_H_

#include <CNSconstants.h>

//////////////////////////////// PaSR TEMPLATE /////////////////////////////////
/**
 * \brief Template for calculating PaSR related subroutines
 *  includes, chemical time scale
 */
template <typename param, typename idx_t>
class PaSR_t {
 
  // input parameters
  amrex::Real C0=1.0;            // Mixing Constant  

  // derived
  //amrex::Real beta =  1.0;


  public:
  PaSR_t() {
     amrex::ParmParse pp("pasr");
      if (!pp.query("mixing_constant", C0)) {
        amrex::Print() << " PaSR: using default mixing constant (1) \n ";   
      }
    
  }

  ~PaSR_t() {}

  /**
   * \brief calculates chemcial time scale for PaSR efficiency model
   * \param[in]  i,j,k        cell index 
   * \param[in]  rY           array of rY after  reaction step
   * \param[in]  Y0           array of Y before reaction step
   * \param[in]  rho          density
   * \param[in]  dt           time step
   * \param[out] tauchem      chemical time scale based on change in species mass fraction during reaction step
   *                          following Domingo school of thought
  */  
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE Real tau_chem(
      const int i, const int j, const int k, const Array4<const Real>& rY, const Array4<const Real>& Y0, 
      const amrex::Real rho, const amrex::Real dt) const
  {

    Real o_dt = 1.0/dt;
    // Calculate laminar chemical source term
    Real omega[NUM_SPECIES];
    for (int n = 0; n < NUM_SPECIES; ++n) {
      omega[n] =  ( rY(i, j, k, n) - rho*Y0(i, j, k, idx_t::QFS + n))* o_dt ;
    }

    // Calculate chemical timescale tau_chem = min(rY/|omega|)
    Real tau = 1e10; // large value if no reaction, can be tuned 
    for (int n = 0; n < NUM_SPECIES; ++n) {
      tau = std::min( rY(i, j, k, n) / std::max( std::abs(omega[n]), 1e-16) , tau);
    }

    return(tau);
  }  
  

////////////////////////////////////////////////////////////////
};



////////////////////////////////////////////////////////////////////////////////
#endif

