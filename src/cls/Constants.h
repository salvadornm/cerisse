#ifndef CONSTANTS_H
#define CONSTANTS_H


using namespace amrex;

namespace universal_constants {
    
    // Speed of light in vacuum (m/s)
    constexpr Real speed_of_light = 299792458.0; // m/s

    // Planck constant (J·s)
    constexpr Real planck_constant = 6.62607015e-34; // J·s

    // Elementary charge (Coulombs)
    constexpr Real elementary_charge = 1.602176634e-19; // C

    // Avogadro constant (mol^-1)
    constexpr Real avogadro_constant = 6.02214076e23; // mol^-1

    // Boltzmann constant (J/K)
    constexpr Real boltzmann_constant = 1.380649e-23; // J·K^-1

    // Universal gas constant (J·mol^-1·K^-1)
    constexpr Real gas_constant = avogadro_constant*boltzmann_constant; // J·mol^-1·K^-1

    // Electron mass (kg)
    constexpr Real electron_mass = 9.10938356e-31; // kg

    // Proton mass (kg)
    constexpr Real proton_mass = 1.67262192369e-27; // kg

    // Neutron mass (kg)
    constexpr Real neutron_mass = 1.67492749804e-27; // kg

    // Pressure units
    constexpr Real pres_bar2si = 1.e5;                  // bar->Pa
    constexpr Real pres_atm2si = 1.01325* pres_bar2si;  // atm->Pa

    // CGS <-> SI (mass)
    constexpr Real mass_cgs2si = 0.001;              // g->kg
    constexpr Real mass_si2cgs = 1.0/mass_cgs2si;    // kg->g
    
    // CGS <-> SI (speed)
    constexpr Real speed_cgs2si = 0.01;              // cm/s->m/s
    constexpr Real speed_si2cgs = 1.0/speed_cgs2si;  // m/s->cm/s
    
    // CGS <-> SI (density) 
    constexpr Real rho_cgs2si = 1000.0;               // g/cm3->kg/m3
    constexpr Real rho_si2cgs = 1.0/rho_cgs2si;       // kg/m3->g/cm3
    
    // CGS <-> SI (pressure)  
    constexpr Real pres_cgs2si = 0.1;                  // dyn/cm2->Pa
    constexpr Real pres_si2cgs = 1.0/pres_cgs2si;      // Pa ->dyn/cm2

    // CGS <-> SI (mass specific energy)  
    constexpr Real specenergy_cgs2si = 1.0e-4;                  // erg/g->J/kg
    constexpr Real specenergy_si2cgs = 1.0/specenergy_cgs2si;   // J/kg->erg/g 

    // CGS <-> SI (volume specific energy rho e)  
    constexpr Real rhoenergy_cgs2si = rho_cgs2si*specenergy_cgs2si;  // erg/cm3->J/m3
    constexpr Real rhoenergy_si2cgs = rho_si2cgs*specenergy_si2cgs;  // J/m3 >erg/cm3 

    // CGS <-> SI (dynamic viscosity)  
    constexpr Real visc_cgs2si = 0.1;                     // Poise->Pa.s  
    constexpr Real visc_si2cgs = 1.0/visc_cgs2si;         // Pa.s->Poise 

    // CGS <-> SI (kinematic viscosity)  
    constexpr Real nuvisc_cgs2si = 1.0e-4;                // Stokes->m2/s
    constexpr Real nuvisc_si2cgs = 1.0/nuvisc_cgs2si;     // m2/s->Stokes 

    // CGS <-> SI (thermal conductivity)  
    constexpr Real cond_cgs2si = 1.0e-6;                  // erg/cm K->W/m K
    constexpr Real cond_si2cgs = 1.0/cond_cgs2si;         // W/m K-> erg/cm K

    // CGS <-> SI (mass diffussivity)  
    constexpr Real diff_cgs2si = 1.0e-4;                  // cm2/s -> m2/s
    constexpr Real diff_si2cgs = 1.0/diff_cgs2si;         // m2/s  -> cm2/s K

    // CGS <-> SI (rho*mass diffussivity) same convesrion as dynamic viscosity  
    constexpr Real rhodiff_cgs2si = 0.1;                   // g/cm s -> kg/m s
    constexpr Real rhodiff_si2cgs = 1.0/rhodiff_cgs2si;    // kg/m s -> g/cm s
     
    // mol <-> kmol
    constexpr Real mol2kmol = 0.001;                 // mol->kmol
    constexpr Real kmol2mol = 1.0/mol2kmol;          // kmol->mol
    
}




#endif // CONSTANTS_H

