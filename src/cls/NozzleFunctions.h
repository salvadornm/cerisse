#ifndef NOZZLEFUNCTIONS_H
#define NOZZLEFUNCTIONS_H

#include <cmath>
#include <array>
#include <numbers>

using namespace amrex;

// Isentropic-flow and nozzle geometry utilities for problem set-up.
// These functions are host-only helpers used to pre-compute initial/boundary
// conditions; they are not called inside GPU kernels.
// Note: the specific gas constant R = 287 J/(kg·K) assumes dry air (γ = 1.4).
//       Pass the appropriate Rgas for other working fluids.
namespace nozzle_functions {

  // Static temperature from stagnation temperature and Mach number (isentropic).
  static inline Real T_isen(Real T0, Real M, Real gamma)
  {
    return T0 / (1.0 + 0.5 * (gamma - 1.0) * M * M);
  }

  // Static pressure from stagnation pressure and Mach number (isentropic).
  static inline Real P_isen(Real P0, Real M, Real gamma)
  {
    const Real factor = 1.0 + 0.5 * (gamma - 1.0) * M * M;
    return P0 * std::pow(factor, -gamma / (gamma - 1.0));
  }

  // Density from stagnation conditions and Mach number (ideal gas, R = 287 J/kg/K).
  static inline Real rho_isen(Real P0, Real T0, Real M, Real gamma,
                               Real Rgas = 287.0)
  {
    const Real T = T_isen(T0, M, gamma);
    const Real P = P_isen(P0, M, gamma);
    return P / (Rgas * T);
  }

  // Area ratio A/A* for isentropic flow at Mach M (choked throat reference).
  static inline Real nozzle_area_ratio(Real M, Real gamma)
  {
    const Real term = (2.0 / (gamma + 1.0)) * (1.0 + 0.5 * (gamma - 1.0) * M * M);
    const Real exponent = (gamma + 1.0) / (2.0 * (gamma - 1.0));
    return std::pow(term, exponent) / M;
  }

  // Stagnation pressure from static pressure, Mach, and γ.
  static inline Real Pstag(Real P, Real M, Real gamma)
  {
    const Real term = 1.0 + 0.5 * (gamma - 1.0) * M * M;
    return P * std::pow(term, gamma / (gamma - 1.0));
  }

  // Stagnation temperature from static temperature, Mach, and γ.
  static inline Real Tstag(Real T, Real M, Real gamma)
  {
    return T * (1.0 + 0.5 * (gamma - 1.0) * M * M);
  }

  // Choked (throat, M=1) pressure from stagnation pressure.
  static inline Real Pchok(Real P0, Real gamma)
  {
    return P0 * std::pow(2.0 / (gamma + 1.0), gamma / (gamma - 1.0));
  }

  // Choked (throat, M=1) temperature from stagnation temperature.
  static inline Real Tchok(Real T0, Real gamma)
  {
    return T0 * (2.0 / (gamma + 1.0));
  }

  // Choked mass-flow rate [kg/s] per unit area (R = 287 J/kg/K for air).
  static inline Real masschok(Real T0, Real P0, Real gamma, Real Rgas = 287.0)
  {
    const Real gp1_o2 = 0.5 * (gamma + 1.0);
    return P0 * std::sqrt(gamma / (Rgas * T0))
               * std::pow(1.0 / gp1_o2, gp1_o2 / (gamma - 1.0));
  }

  // Static pressure at Mach M in a nozzle with stagnation pressure P0.
  static inline Real Pnozz(Real Mach, Real P0, Real gamma)
  {
    if (Mach < 0.0) return 0.0;
    const Real factor = 1.0 + 0.5 * (gamma - 1.0) * Mach * Mach;
    return P0 * std::pow(factor, -gamma / (gamma - 1.0));
  }

  // Nozzle-exit location on a sphere-cone geometry.
  // Returns {x, y, z} in the body-fixed frame given:
  //   factor_y        : fractional radial position on cone surface
  //   theta_nz_deg    : azimuthal nozzle angle [degrees]
  //   radius_cone     : base radius of the cone
  //   theta_cone_deg  : half-angle of the cone [degrees]
  //   R_nozzle_exit   : non-dimensional nozzle-exit radius (relative to radius_cone)
  static inline std::array<Real, 3> compute_srp_xyz(
      Real factor_y, Real theta_nz_deg,
      Real radius_cone, Real theta_cone_deg, Real R_nozzle_exit)
  {
    const Real pi = std::numbers::pi;
    const Real theta_cone = theta_cone_deg * (pi / 180.0);
    const Real theta_nz   = theta_nz_deg   * (pi / 180.0);

    const Real h_cone = radius_cone / std::tan(theta_cone);
    const Real x0     = h_cone * R_nozzle_exit / radius_cone - h_cone;
    const Real y0     = factor_y * radius_cone;

    const Real x = x0 + y0 / std::tan(theta_cone);
    const Real y = y0 * std::cos(theta_nz);
    const Real z = y0 * std::sin(theta_nz);

    return {x, y, z};
  }

} // namespace nozzle_functions

#endif // NOZZLEFUNCTIONS_H
