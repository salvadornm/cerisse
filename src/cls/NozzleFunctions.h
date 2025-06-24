#ifndef NOZZLEFUNCTIONS_H
#define NOZZLEFUNCTIONS_H

#include <cmath>

using namespace amrex;


namespace nozzle_functions {
  // Set of useful functions to set-up nozzle cases
  // These are not used in any calculations within Cerisse 
  // they  are only to help 
  // (aditionally use nozzle_calculator in cerisse/tools to compute nozzle parameters)

  static constexpr Real T_isen(const Real T0, const Real M, const Real gamma)
  {
    return T0 / (1.0 + 0.5 * (gamma - 1.0) * M * M);
  }

  static constexpr Real P_isen(const Real P0, const Real M, const Real gamma)
  {
    const Real factor = 1.0 + 0.5 * (gamma - 1.0) * M * M;
    return P0 * std::pow(factor, -gamma / (gamma - 1.0));
  }

  static constexpr Real rho_isen(const Real P0, const Real T0, const Real M, const Real gamma)
  {
    const Real R = 287; // for air
    const Real T = T_isen(T0, M, gamma);
    const Real P = P_isen(P0, M, gamma);
    return P / (R * T);
  }

  static constexpr Real nozzle_area_ratio(const Real M, const Real gamma)
{
    Real term1 = 1.0 / M;
    Real term2 = (2.0 / (gamma + 1.0)) * (1.0 + 0.5 * (gamma - 1.0) * M * M);
    Real exponent = (gamma + 1.0) / (2.0 * (gamma - 1.0));
    return term1 * std::pow(term2, exponent);
}

  // \brief computes stagnation pressure,function P,M,gamma
  // \param P:     pressure
  // \param M:     Mach number
  // \param gamma: gamma (assumes ideal gas)  
  static const Real Pstag(const Real P, const Real M, const Real gamma)
  {         
    const Real term = 1.0 + 0.5 * (gamma - 1.0) * M * M;
    const Real exponent = gamma / (gamma - 1.0);
    return P * std::pow(term, exponent);
  }
  // \brief computes stagnation temperature,function T,M,gamma
  // \param T:     Temperature
  // \param M:     Mach number
  // \param gamma: gamma (assumes ideal gas)  
  static const Real Tstag(const Real T, const Real M, const Real gamma)
  {         
    Real T0 = T*( 1 + 0.5*(gamma-1)*M*M) ;
    return(T0);
  }
  // \brief Computes the choked pressure (throat pressure) for isentropic flow
  // Assumes ideal gas and Mach = 1
  // \param P0:    stagnation pressure  
  // \param gamma: gamma  
  static const Real Pchok(const Real P0, const Real gamma)
  {
    const Real exponent = gamma / (gamma - 1.0);
    return P0 * std::pow(2.0 / (gamma + 1.0), exponent);
  }
  // Computes the choked temperature (throat temperature) for isentropic flow
  // Assumes ideal gas and Mach = 1
  // \param P0:    stagnation pressure  
  // \param gamma: gamma  
  static const Real Tchok(const Real T0, const Real gamma)
  {
    return T0 * (2.0 / (gamma + 1.0));
  }  
  // \brief computes choked mass flow rate [kg/s]
  // \param T0:     stagnation Temperature
  // \param P0:     stagnation pressure
  // \param gamma: gamma (assumes ideal gas)
  static const Real masschok(const Real T0, const Real P0, const Real gamma)
  {      
    const Real R = 287; //air
    const Real gamma_p1_o2 = 0.5*(gamma+1.0);
    const Real m = P0*sqrt(gamma/(R*T0))*std::pow(1.0/gamma_p1_o2,gamma_p1_o2/(gamma-1));
    return(m);
  }
  // \brief computes Presure as function Ma,P0 and gam
  // \param P0:     stagnation pressure
  static const Real Pnozz(const Real Mach, const Real P0, const Real gamma) {
    if (Mach < 0.0) {return 0.0;}
    const Real factor = 1.0 + 0.5*(gamma - 1.0) * Mach * Mach;
    return P0 * std::pow(factor, -gamma / (gamma - 1.0));
  }

}

//// calculate nozzle location on a sphere cone
static constexpr std::array<Real, 3> compute_srp_xyz(Real factor_y, Real theta_nz_deg, Real radius_cone, Real theta_cone_deg, Real R_nozzle_exit)
{
  Real theta_cone_rad = theta_cone_deg * (std::numbers::pi / 180.0);
  Real theta_nz_rad   = theta_nz_deg * (std::numbers::pi / 180.0);

  Real h_cone = radius_cone/std::tan(theta_cone_rad);
  Real x0     = h_cone*R_nozzle_exit/radius_cone - h_cone;
  Real y0     = factor_y * radius_cone;

  Real x      = x0 + y0 / std::tan(theta_cone_rad);
  Real y      = y0 * std::cos(theta_nz_rad);
  Real z      = y0 * std::sin(theta_nz_rad);

  return {x, y, z};
} 

#endif // NOZZLEFUNCTIONS_H

