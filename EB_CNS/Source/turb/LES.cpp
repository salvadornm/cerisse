#include "LES.H"

using namespace amrex;

void LESModel::create(std::string selector, LESModel*& model)
{
  LESModel** tmp = (LESModel**)The_Arena()->alloc(sizeof(LESModel*));
  if (selector == Smagorinsky::identifier()) {
    amrex::single_task([=] AMREX_GPU_DEVICE() { *tmp = new Smagorinsky(); });
  } else if (selector == WALE::identifier()) {
    amrex::single_task([=] AMREX_GPU_DEVICE() { *tmp = new WALE(); });
  } else {
    amrex::Abort("Cannot find " + selector + " in " + base_identifier());
  }
  amrex::Gpu::copy(Gpu::deviceToHost, tmp, tmp + 1, &model);
  The_Arena()->free(tmp);
}

void LESModel::destroy(LESModel*& model) {
  auto local_model = model;
  amrex::single_task([=] AMREX_GPU_DEVICE() { delete local_model; });
}

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

/**
 * @brief Calculate turbulent viscosity using the Smagorinsky model.
 *        Note that velocity derivatives are evaluated at cell centre (cc).
 *
 * @param i,j,k     x,y,z index
 * @param q         primitive variables array
 * @param dxinv     cell size, used for calculating velocity derivatives
 * @param deltabar  filter width
 * @param Cs        Smagorinsky constant
 * @param mu_T[out] output turbulent viscosity
 */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void Smagorinsky::mu_T_cc(
  const int i, const int j, const int k, const Array4<const Real>& q,
  const GpuArray<Real, AMREX_SPACEDIM>& dxinv, const Real deltabar, const Real Cs,
  Real& mu_T)
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

  Real Sijmag = 0;
  for (int m = 0; m < AMREX_SPACEDIM; m++) {
    for (int n = 0; n < AMREX_SPACEDIM; n++) {
      Sijmag += (0.5 * (dUdx[m][n] + dUdx[n][m])) *
                (0.5 * (dUdx[m][n] + dUdx[n][m])); // Sij*Sij
    }
  }
  Sijmag = std::sqrt(2.0 * Sijmag);

  mu_T = q(i, j, k, QRHO) * Cs * Cs * deltabar * deltabar * Sijmag;
}

/**
 * @brief Calculate turbulent viscosity using the Smagorinsky model.
 *        Note that velocity derivatives are evaluated at face centre (fc).
 *
 * @param dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz velocity derivatives at face
 * @param rho       density at face
 * @param deltabar  filter width
 * @param Cs        Smagorinsky constant
 * @param mu_T[out] output turbulent viscosity divided by Cs^2
 */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void Smagorinsky::mu_T_fc(
  const Real dudx, const Real dudy, const Real dudz, const Real dvdx,
  const Real dvdy, const Real dvdz, const Real dwdx, const Real dwdy,
  const Real dwdz, const Real rho, const Real deltabar, const Real Cs, Real& mu_T)
{
  Real dUdx[3][3] = {{0.0}};
  dUdx[0][0] = dudx;
  dUdx[1][0] = dvdx;
  dUdx[0][1] = dudy;
  dUdx[1][1] = dvdy;
#if AMREX_SPACEDIM == 3
  dUdx[2][0] = dwdx;
  dUdx[2][1] = dwdy;
  dUdx[0][2] = dudz;
  dUdx[1][2] = dvdz;
  dUdx[2][2] = dwdz;
#endif

  Real Sijmag = 0;
  for (int m = 0; m < AMREX_SPACEDIM; m++) {
    for (int n = 0; n < AMREX_SPACEDIM; n++) {
      Sijmag += (0.5 * (dUdx[m][n] + dUdx[n][m])) *
                (0.5 * (dUdx[m][n] + dUdx[n][m])); // Sij*Sij
    }
  }
  Sijmag = std::sqrt(2.0 * Sijmag);

  mu_T = rho * Cs * Cs * deltabar * deltabar * Sijmag;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
WALE::mu_T_cc(const int i, const int j, const int k, const Array4<const Real>& q,
              const GpuArray<Real, AMREX_SPACEDIM>& dxinv, const Real deltabar,
              const Real Cs, Real& mu_T)
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

  const Real divu = AMREX_D_TERM(dUdx[0][0], +dUdx[1][1], +dUdx[2][2]);
  Real dUdx2[3][3] = {{0.0}};
  for (int m = 0; m < AMREX_SPACEDIM; m++) {
    for (int n = 0; n < AMREX_SPACEDIM; n++) {
      dUdx2[m][n] =
        dUdx[m][0] * dUdx[0][n] + dUdx[m][1] * dUdx[1][n] + dUdx[m][2] * dUdx[2][n];
    }
  }
  // assert(divu * divu == dUdx2[0][0] + dUdx2[1][1] + dUdx2[2][2]); // tested true

  Real SijSij = 0.0;
  Real DijDij = 0.0;
  for (int m = 0; m < AMREX_SPACEDIM; m++) {
    for (int n = 0; n < AMREX_SPACEDIM; n++) {
      SijSij += (0.5 * (dUdx[m][n] + dUdx[n][m])) *
                (0.5 * (dUdx[m][n] + dUdx[n][m])); // Sij*Sij
      DijDij +=
        (0.5 * (dUdx2[m][n] + dUdx2[n][m]) - Real(m == n) * divu * divu / 3) *
        (0.5 * (dUdx2[m][n] + dUdx2[n][m]) -
         Real(m == n) * divu * divu / 3); // Dij*Dij
    }
  }

  Real Cw = std::sqrt(10.6) * Cs;
  mu_T = q(i, j, k, QRHO) * Cw * Cw * deltabar * deltabar * std::pow(DijDij, 1.5) /
         (std::pow(SijSij, 2.5) + std::pow(DijDij, 1.25) +
          std::numeric_limits<Real>::denorm_min());
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void WALE::mu_T_fc(
  const Real dudx, const Real dudy, const Real dudz, const Real dvdx,
  const Real dvdy, const Real dvdz, const Real dwdx, const Real dwdy,
  const Real dwdz, const Real rho, const Real deltabar, const Real Cs, Real& mu_T)
{
  Real dUdx[3][3] = {{0.0}};
  dUdx[0][0] = dudx;
  dUdx[1][0] = dvdx;
  dUdx[0][1] = dudy;
  dUdx[1][1] = dvdy;
#if AMREX_SPACEDIM == 3
  dUdx[2][0] = dwdx;
  dUdx[2][1] = dwdy;
  dUdx[0][2] = dudz;
  dUdx[1][2] = dvdz;
  dUdx[2][2] = dwdz;
#endif

  const Real divu = AMREX_D_TERM(dUdx[0][0], +dUdx[1][1], +dUdx[2][2]);

  Real Sijmag = 0;
  Real Dijmag = 0;
  for (int m = 0; m < AMREX_SPACEDIM; m++) {
    for (int n = 0; n < AMREX_SPACEDIM; n++) {
      Sijmag += (0.5 * (dUdx[m][n] + dUdx[n][m])) *
                (0.5 * (dUdx[m][n] + dUdx[n][m])); // Sij*Sij
      Dijmag +=
        (0.5 * (dUdx[m][n] * dUdx[m][n] + dUdx[n][m] * dUdx[n][m]) -
         Real(m == n) * divu * divu / 3) *
        (0.5 * (dUdx[m][n] * dUdx[m][n] + dUdx[n][m] * dUdx[n][m]) -
         Real(m == n) * divu * divu / 3); // Dij*Dij // TODO: Is this correct?
    }
  }

  Real Cw = std::sqrt(10.6) * Cs;
  mu_T = rho * Cw * Cw * deltabar * deltabar * std::pow(Dijmag, 1.5) /
         (std::pow(Sijmag, 2.5) + std::pow(Dijmag, 1.25));
}