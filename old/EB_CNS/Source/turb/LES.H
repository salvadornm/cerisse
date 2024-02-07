#ifndef CNS_LES_H_
#define CNS_LES_H_

#include <AMReX_Array4.H>
#include <AMReX_REAL.H>
#include <Factory.H>

#include "index_macros.H"

using namespace amrex;

/**
 * @brief Base class to register eddy viscosity type LES models.
 */
class LESModel : public pele::physics::Factory<LESModel>
{
public:
  static const std::string base_identifier() { return "LESModel"; }

  ~LESModel() override = default;

  virtual void mu_T_cc(const int i, const int j, const int k,
                       const Array4<const Real>& q,
                       const GpuArray<Real, AMREX_SPACEDIM>& dxinv,
                       const Real deltabar, const Real Cs, Real& mu_T) = 0;

  virtual void mu_T_fc(const Real dudx, const Real dudy, const Real dudz,
                       const Real dvdx, const Real dvdy, const Real dvdz,
                       const Real dwdx, const Real dwdy, const Real dwdz,
                       const Real rho, const Real deltabar, const Real Cs,
                       Real& mu_T) = 0;

  virtual void xi_T_cc(const int i, const int j, const int k,
                       const Array4<const Real>& q,
                       const GpuArray<Real, AMREX_SPACEDIM>& dxinv,
                       const Real deltabar, const Real CI, Real& xi_T);
};

/**
 * @brief Calculate the turbulent bulk viscosity.
 *        Note that velocity derivatives are evaluated at cell centre (cc).
 *
 * @param i,j,k     x,y,z index
 * @param q         primitive variables array
 * @param dxinv     cell size, used for calculating velocity derivatives
 * @param deltabar  filter width
 * @param CI        Yoshizawa constant
 * @param xi_T[out] output turbulent bulk viscosity
 */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
LESModel::xi_T_cc(const int i, const int j, const int k, const Array4<const Real>& q,
                  const GpuArray<Real, AMREX_SPACEDIM>& dxinv, const Real deltabar,
                  const Real CI, Real& xi_T)
{
  // Calculate derivatives at cell centers, second order central difference
  Real dUdx[3][3] = {{0.0}};
  dUdx[0][0] = 0.5 * dxinv[0] * (q(i + 1, j, k, QU) - q(i - 1, j, k, QU)); // dudx
  dUdx[1][0] = 0.5 * dxinv[0] * (q(i + 1, j, k, QV) - q(i - 1, j, k, QV)); // dvdx
  dUdx[0][1] = 0.5 * dxinv[1] * (q(i, j + 1, k, QU) - q(i, j - 1, k, QU)); // dudy
  dUdx[1][1] = 0.5 * dxinv[1] * (q(i, j + 1, k, QV) - q(i, j - 1, k, QV)); // dvdy
#if AMREX_SPACEDIM == 3
  dUdx[2][0] = 0.5 * dxinv[0] * (q(i + 1, j, k, QW) - q(i - 1, j, k, QW)); // dwdx
  dUdx[2][1] = 0.5 * dxinv[1] * (q(i, j + 1, k, QW) - q(i, j - 1, k, QW)); // dwdy
  dUdx[0][2] = 0.5 * dxinv[2] * (q(i, j, k + 1, QU) - q(i, j, k - 1, QU)); // dudz
  dUdx[1][2] = 0.5 * dxinv[2] * (q(i, j, k + 1, QV) - q(i, j, k - 1, QV)); // dvdz
  dUdx[2][2] = 0.5 * dxinv[2] * (q(i, j, k + 1, QW) - q(i, j, k - 1, QW)); // dwdz
#endif

  Real Sijmag2 = 0.0, div = 0.0;
  for (int m = 0; m < AMREX_SPACEDIM; m++) {
    for (int n = 0; n < AMREX_SPACEDIM; n++) {
      Sijmag2 += (0.5 * (dUdx[m][n] + dUdx[n][m])) *
                 (0.5 * (dUdx[m][n] + dUdx[n][m])); // Sij*Sij
    }
    div += dUdx[m][m]; // Skk
  }
  Sijmag2 *= 2.0;

  if (std::abs(div) > std::numeric_limits<Real>::epsilon()) {
    xi_T = 1.0 / 3.0 * q(i, j, k, QRHO) * CI * deltabar * deltabar * Sijmag2 / div;
  } else {
    xi_T = 0.0;
  }
}

//
// Declare LES models

class Smagorinsky : public LESModel::Register<Smagorinsky>
{
public:
  static const std::string identifier() { return "Smagorinsky"; }

  void mu_T_cc(const int i, const int j, const int k, const Array4<const Real>& q,
               const GpuArray<Real, AMREX_SPACEDIM>& dxinv, const Real deltabar,
               const Real Cs, Real& mu_T) override;

  void mu_T_fc(const Real dudx, const Real dudy, const Real dudz, const Real dvdx,
               const Real dvdy, const Real dvdz, const Real dwdx, const Real dwdy,
               const Real dwdz, const Real rho, const Real deltabar, const Real Cs,
               Real& mu_T) override;
};

// Nicoud & Ducros (1999) https://doi.org/10.1023/A:1009995426001
class WALE : public LESModel::Register<WALE>
{
public:
  static const std::string identifier() { return "WALE"; }

  void mu_T_cc(const int i, const int j, const int k, const Array4<const Real>& q,
               const GpuArray<Real, AMREX_SPACEDIM>& dxinv, const Real deltabar,
               const Real Cs, Real& mu_T) override;

  void mu_T_fc(const Real dudx, const Real dudy, const Real dudz, const Real dvdx,
               const Real dvdy, const Real dvdz, const Real dwdx, const Real dwdy,
               const Real dwdz, const Real rho, const Real deltabar, const Real Cs,
               Real& mu_T) override;
};

#endif