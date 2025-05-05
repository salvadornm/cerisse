#ifndef COMBUSTIONFUNCTIONS_H
#define COMBUSTIONFUNCTIONS_H


using namespace amrex;


namespace combustion_functions {
  // Set of useful functions to set-up combustion cases
  // These are not use for any calculations within Cerisse 
  // they  are only to help

  constexpr Real molarMassH2  = 2.016;  // [kg/kmol]
  constexpr Real molarMassAir = 28.97;
  constexpr Real molarMassO2  = 31.999;
  constexpr Real molarMassN2  = 28.0134;    
  constexpr Real XO2air = 0.21;
  constexpr Real XN2air = 0.79;

  struct Moles_Mixture
  {
    Real H2;
    Real N2;
    Real O2;
  };
  struct Y_Mixture
  {
    Real H2;
    Real N2;
    Real O2;
  };

  inline Moles_Mixture compute_MixMol_fromEquivalenceRatio_H2Air(Real phi)
  {      
    const Real AFR_stoic = 132.476/4.032;
    Real AFR     = AFR_stoic/phi;
    Real air_mol = AFR*molarMassAir/molarMassH2;
    return{1.0,air_mol*XO2air,air_mol*XN2air};
  }

  inline Y_Mixture mol2Y_H2Air(Moles_Mixture mix)
  {
    Real MMW = mix.H2*molarMassH2 +mix.O2*molarMassO2 + mix.N2*molarMassN2;  
    return{mix.H2*molarMassH2/MMW,mix.N2*molarMassN2/MMW,mix.O2*molarMassO2/MMW};
  } 

}

#endif // COMBUSTIONFUNCTIONS_H

