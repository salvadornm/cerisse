#ifndef NOZZLEFUNCTIONS_JIAYE_H
#define NOZZLEFUNCTIONS_JIAYE_H

#include <cmath>
#include <numbers>

using namespace amrex;


namespace nozzle_functions {

class Nozzle {

public:
  const int Nozzle_Type;
  bool Nozzle_Switch;
  const std::array<Real, AMREX_SPACEDIM> Exit_Centre;
  const std::array<Real, AMREX_SPACEDIM> Orientation;
  const Real Mach_Exit;
  Real P_Stag;
  Real T_Stag;
  const Real R_Exit;
  const Real R_Plenum;
  const Real Divergence_Deg;
  const Real Convergence_Deg;
  const Real Throat_Length;
  const Real Plenum_Length;
  const Real gamma;
  const Real Rgas;

  constexpr Nozzle(
            int  Nozzle_Type_,
            int  Nozzle_Switch_,
            std::array<Real,AMREX_SPACEDIM> Exit_Centre_,
            std::array<Real,AMREX_SPACEDIM> Orientation_,
            Real Mach_Exit_,
            Real P_Stag_,
            Real T_Stag_,
            Real R_Exit_,
            Real R_Plenum_         = 0.0,
            Real Divergence_Deg_   = 15.0,
            Real Convergence_Deg_  = 20.0,
            Real Throat_Length_    = 0.0,
            Real Plenum_Length_    = 0.0,
            Real gamma_            = 1.4,
            Real Rgas_             = 287.0512)
  :
  Nozzle_Type(Nozzle_Type_)           ,   Nozzle_Switch(Nozzle_Switch_)       ,
  Exit_Centre(Exit_Centre_)           ,   Orientation(Orientation_)           ,
  Mach_Exit(Mach_Exit_)               ,        
  P_Stag(P_Stag_)                     ,   T_Stag(T_Stag_)                     , 
  R_Exit(R_Exit_)                     ,   R_Plenum(R_Plenum_)                 ,
  Divergence_Deg(Divergence_Deg_)     ,   Convergence_Deg(Convergence_Deg_)   , 
  Throat_Length(Throat_Length_)       ,   Plenum_Length(Plenum_Length_)       ,
  gamma(gamma_)                       ,   Rgas(Rgas_)                        {}

  constexpr Real Aera_Ratio() const {
    const Real term1    = 1.0 / Mach_Exit;
    const Real term2    = (2.0 / (gamma + 1.0)) * (1.0 + 0.5 * (gamma - 1.0) * Mach_Exit * Mach_Exit);
    const Real exponent = (gamma + 1.0) / (2.0 * (gamma - 1.0));
    return term1 * std::pow(term2, exponent); }

  constexpr Real R_Throat () const {
    return R_Exit / std::sqrt(Aera_Ratio()); }

  constexpr Real A_Exit() const {
    return std::numbers::pi * R_Exit * R_Exit; }

  constexpr Real A_Throat() const {
    return std::numbers::pi * R_Throat () * R_Throat (); }

  constexpr Real A_Plenum() const {
    return std::numbers::pi * R_Plenum * R_Plenum; } 

  constexpr Real Divergence_Length() const {
  return (R_Exit - R_Throat ()) / std::tan(Divergence_Deg * std::numbers::pi/180.0); }

  constexpr Real Convergence_Length() const {
    return (std::max(R_Plenum, R_Throat()) - R_Throat ()) / std::tan(Convergence_Deg * std::numbers::pi/180.0); }

  constexpr Real Nozzle_Length() const {
    switch (Nozzle_Type) {
      case 0:
        return 0.0;
      case 1:
        return Divergence_Length();
      case 2:
        return Divergence_Length() + Throat_Length + Convergence_Length();
      default:
        return Divergence_Length(); } }


  constexpr Real P_Exit() const {
    const Real factor = 1.0 + 0.5 * (gamma - 1.0) * Mach_Exit * Mach_Exit;
    return P_Stag * std::pow(factor, -gamma / (gamma - 1.0)); }

  constexpr Real T_Exit() const {
    return T_Stag / (1.0 + 0.5 * (gamma - 1.0) * Mach_Exit * Mach_Exit); }

  constexpr Real Rho_Exit() const {
    return P_Exit() / (Rgas * T_Exit()); }

  constexpr Real V_Exit() const {
    return Mach_Exit * std::sqrt(gamma * Rgas * T_Exit()); }


  constexpr Real P_Throat() const {
    const Real exponent = gamma / (gamma - 1.0);
    return P_Stag * std::pow(2.0 / (gamma + 1.0), exponent); }

  constexpr Real T_Throat() const {
    return T_Stag * (2.0 / (gamma + 1.0)); }

  constexpr Real Rho_Throat() const {
    return P_Throat() / (Rgas * T_Throat()); }

  constexpr Real V_Throat() const { 
    return std::sqrt(gamma * Rgas * T_Throat()); }

  constexpr Real Rho_Stag() const {
    return P_Stag / (Rgas * T_Stag); }
    

  constexpr Real Mass_flow() const { 
    return std::sqrt(gamma * Rgas * T_Throat()); }

  constexpr Real Thrust() const { 
    return std::sqrt(gamma * Rgas * T_Throat()); }


  constexpr Real R_BC_Applied() const {
    switch (Nozzle_Type) {
      case 0:
        return R_Exit;
      case 1:
        return R_Throat();
      case 2:
        return R_Plenum;
      default:
        return R_Throat(); } }

  constexpr Real P_BC_Applied() const {
    switch (Nozzle_Type) {
      case 0:
        return P_Exit();
      case 1:
        return P_Throat();
      case 2:
        return P_Stag;
      default:
        return P_Throat(); } }

  constexpr Real T_BC_Applied() const {
    switch (Nozzle_Type) {
      case 0:
        return T_Exit();
      case 1:
        return T_Throat();
      case 2:
        return T_Stag;
      default:
        return T_Throat(); } } 
        
  constexpr Real Rho_BC_Applied() const {
    switch (Nozzle_Type) {
      case 0:
        return Rho_Exit();
      case 1:
        return Rho_Throat();
      case 2:
        return Rho_Stag();
      default:
        return Rho_Throat(); } }   
        
  constexpr Real V_BC_Applied() const {
    switch (Nozzle_Type) {
      case 0:
        return V_Exit();
      case 1:
        return V_Throat();
      case 2:
        return 0.0;
      default:
        return V_Throat(); } } 


  constexpr Real Orientation_Norm() const {
    switch (AMREX_SPACEDIM) {
      case 2:
        return std::sqrt(Orientation[0]*Orientation[0] + Orientation[1]*Orientation[1]);
      case 3:
        return std::sqrt(Orientation[0]*Orientation[0] + Orientation[1]*Orientation[1] + Orientation[2]*Orientation[2]);
      default:
        return std::sqrt(Orientation[0]*Orientation[0] + Orientation[1]*Orientation[1] + Orientation[2]*Orientation[2]); }
  }

  constexpr std::array<Real, AMREX_SPACEDIM> BC_Applied_Centre() const{
    Real x = Exit_Centre[0] - Nozzle_Length() * Orientation[0] / Orientation_Norm();
    Real y = Exit_Centre[1] - Nozzle_Length() * Orientation[1] / Orientation_Norm();
    switch (AMREX_SPACEDIM) {
      case 2:
        return {x, y};
      case 3:  
        return {x, y, Exit_Centre[2] - Nozzle_Length() * Orientation[2] / Orientation_Norm()};
      default: 
        return {x, y, Exit_Centre[2] - Nozzle_Length() * Orientation[2] / Orientation_Norm()}; }
  } 

  constexpr std::array<Real, AMREX_SPACEDIM> Nozzle_Centre() const{
    Real x = Exit_Centre[0] - 0.5*Nozzle_Length() * Orientation[0] / Orientation_Norm();
    Real y = Exit_Centre[1] - 0.5*Nozzle_Length() * Orientation[1] / Orientation_Norm();
    switch (AMREX_SPACEDIM) {
      case 2:
        return {x, y};
      case 3:  
        return {x, y, Exit_Centre[2] - 0.5*Nozzle_Length() * Orientation[2] / Orientation_Norm()};
      default: 
        return {x, y, Exit_Centre[2] - 0.5*Nozzle_Length() * Orientation[2] / Orientation_Norm()}; }
  } 

  constexpr Real Nozzle_Vol() const {
    Real Width = std::max(R_Exit, R_Plenum);
    return Width * Width + 0.5*Nozzle_Length() * 0.5*Nozzle_Length();
  }


  constexpr void Nozzle_ON()  { Nozzle_Switch = true ; }
  constexpr void Nozzle_OFF() { Nozzle_Switch = false; }

  constexpr void Set_P_Stag(Real P_Stag_) { P_Stag = P_Stag_; }

};

  // Set of useful functions to set-up nozzle cases
  // These are not use for any calculations within Cerisse 
  // they  are only to help (required constexpr)

  static constexpr Real T_isen(const Real T0, const Real M, const Real gammama)
  {
    return T0 / (1.0 + 0.5 * (gammama - 1.0) * M * M);
  }

  static constexpr Real P_isen(const Real P0, const Real M, const Real gammama)
  {
    const Real factor = 1.0 + 0.5 * (gammama - 1.0) * M * M;
    return P0 * std::pow(factor, -gammama / (gammama - 1.0));
  }

  static constexpr Real rho_isen(const Real P0, const Real T0, const Real M, const Real gammama)
  {
    const Real R = 287; // for air
    const Real T = T_isen(T0, M, gammama);
    const Real P = P_isen(P0, M, gammama);
    return P / (R * T);
  }

  static constexpr Real nozzle_area_ratio(const Real M, const Real gammama)
  {
    Real term1 = 1.0 / M;
    Real term2 = (2.0 / (gammama + 1.0)) * (1.0 + 0.5 * (gammama - 1.0) * M * M);
    Real exponent = (gammama + 1.0) / (2.0 * (gammama - 1.0));
    return term1 * std::pow(term2, exponent);
  }

  // \brief computes stagnation pressure,function P,M,gammama
  // \param P:     pressure
  // \param M:     Mach number
  // \param gammama: gammama (assumes ideal gas)  
  static constexpr Real Pstag(const Real P, const Real M, const Real gammama)
  {         
    const Real term = 1.0 + 0.5 * (gammama - 1.0) * M * M;
    const Real exponent = gammama / (gammama - 1.0);
    return P * std::pow(term, exponent);
  }
  // \brief computes stagnation temperature,function T,M,gammama
  // \param T:     Temperature
  // \param M:     Mach number
  // \param gammama: gammama (assumes ideal gas)  
  static constexpr Real Tstag(const Real T, const Real M, const Real gammama)
  {         
    Real T0 = T*( 1 + 0.5*(gammama-1)*M*M) ;
    return(T0);
  }
  // \brief Computes the choked pressure (Throat pressure) for isentropic flow
  // Assumes ideal gas and Mach = 1
  // \param P0:    stagnation pressure  
  // \param gammama: gammama  
  static constexpr Real Pchok(const Real P0, const Real gammama)
  {
    const Real exponent = gammama / (gammama - 1.0);
    return P0 * std::pow(2.0 / (gammama + 1.0), exponent);
  }
  // Computes the choked temperature (Throat temperature) for isentropic flow
  // Assumes ideal gas and Mach = 1
  // \param P0:    stagnation pressure  
  // \param gammama: gammama  
  static constexpr Real Tchok(const Real T0, const Real gammama)
  {
    return T0 * (2.0 / (gammama + 1.0));
  }  
  // \brief computes choked mass flow rate [kg/s]
  // \param T0:     stagnation Temperature
  // \param P0:     stagnation pressure
  // \param gammama: gammama (assumes ideal gas)
  static constexpr Real masschok(const Real T0, const Real P0, const Real gammama)
  {      
    const Real R = 287; //air
    const Real gammama_p1_o2 = 0.5*(gammama+1.0);
    const Real m = P0*sqrt(gammama/(R*T0))*std::pow(1.0/gammama_p1_o2,gammama_p1_o2/(gammama-1));
    return(m);
  }
  // \brief computes Presure as function Ma,P0 and gamma
  // \param P0:     stagnation pressure
  static constexpr Real Pnozz(const Real Mach, const Real P0, const Real gammama) {
    if (Mach < 0.0) {return 0.0;}
    const Real factor = 1.0 + 0.5*(gammama - 1.0) * Mach * Mach;
    return P0 * std::pow(factor, -gammama / (gammama - 1.0));
  }

}

//// calculate nozzle location on a sphere cone
static constexpr std::array<Real, 3> compute_srp_xyz(Real factor_y, Real theta_nz_deg, Real radius_cone, Real theta_cone_deg, Real R_nozzle_Exit)
{
  Real theta_cone_rad = theta_cone_deg * (std::numbers::pi / 180.0);
  Real theta_nz_rad   = theta_nz_deg * (std::numbers::pi / 180.0);

  Real h_cone = radius_cone/std::tan(theta_cone_rad);
  Real x0     = h_cone*R_nozzle_Exit/radius_cone - h_cone;
  Real y0     = factor_y * radius_cone;

  Real x      = x0 + y0 / std::tan(theta_cone_rad);
  Real y      = y0 * std::cos(theta_nz_rad);
  Real z      = y0 * std::sin(theta_nz_rad);

  return {x, y, z};
} 

#endif // NOZZLEFUNCTIONSJIAYE_H

