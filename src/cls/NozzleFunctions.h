#ifndef NOZZLEFUNCTIONS_H
#define NOZZLEFUNCTIONS_H

#include <cmath>

using namespace amrex;


namespace nozzle_functions {
  // Set of useful functions to set-up nozzle cases
  // These are not use for any calculations within Cerisse 
  // they  are only to help (required constexpr)


  // \brief computes stagnation pressure,function P,M,gamma
  // \param P:     pressure
  // \param M:     Mach number
  // \param gamma: gamma (assumes ideal gas)  
  static constexpr Real Pstag(const Real P, const Real M, const Real gamma)
  {         
    const Real term = 1.0 + 0.5 * (gamma - 1.0) * M * M;
    const Real exponent = gamma / (gamma - 1.0);
    return P * std::pow(term, exponent);
  }

  // \brief computes stagnation temperature,function T,M,gamma
  // \param T:     Temperature
  // \param M:     Mach number
  // \param gamma: gamma (assumes ideal gas)  
  static constexpr Real Tstag(const Real T, const Real M, const Real gamma)
  {      
   
    Real T0 = T*( 1 + 0.5*(gamma-1)*M*M) ;

    return(T0);
  }

  // \brief Computes the choked pressure (throat pressure) for isentropic flow
  // Assumes ideal gas and Mach = 1
  // \param P0:    staganation pressure  
  // \param gamma: gamma  
  static constexpr Real Pchok(const Real P0, const Real gamma)
  {
    const Real exponent = gamma / (gamma - 1.0);
    return P0 * std::pow(2.0 / (gamma + 1.0), exponent);
  }
  // Computes the choked temperature (throat temperature) for isentropic flow
  // Assumes ideal gas and Mach = 1
  // \param P0:    staganation pressure  
  // \param gamma: gamma  
  static constexpr Real Tchok(const Real T0, const Real gamma)
  {
    return T0 * (2.0 / (gamma + 1.0));
  }  
  // \brief computes choked mass flow rate [kg/s]
  // \param T0:     stganation Temperature
  // \param P0:     stganation pressure
  // \param gamma: gamma (assumes ideal gas)
  static constexpr Real masschok(const Real T0, const Real P0, const Real gamma)
  {      
    Real R = 287; //air
    Real gamma_p1_o2 = 0.5*(gamma+1.0);
    Real m = P0*sqrt(gamma/(R*T0))*std::pow(1.0/gamma_p1_o2,gamma_p1_o2/(gamma-1));

    return(m);
  }

}

#endif // NOZZLEFUNCTIONS_H

