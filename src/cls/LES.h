#ifndef LES_H_
#define LES_H_

////////////////////////////////TRANSPORT/////////////////////////////////
template <typename param>
class LES_t {
  private:
  Real Prsgs = param::Prsgs;
  Real Scsgs = param::Scsgs;
  
  public:
  

  // calculate derivatives and store in dUdx[n][m] gradvel (QV pass class)
  // smag   (dUdx,Delta)
  // wale   (dUdx,Delta)
  // tausgs (dUdx,Delta)

  AMREX_GPU_DEVICE AMREX_FORCE_INLINE void gradvel(
  const int i, const int j, const int k, const Array4<const Real>& q,
  const GpuArray<Real, AMREX_SPACEDIM>& dxinv, const GpuArray<Real, AMREX_SPACEDIM>& dUdx)
  {

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
 
  }

  AMREX_GPU_DEVICE AMREX_FORCE_INLINE Real Smagorinsky::viscsgs(
  const GpuArray<Real, AMREX_SPACEDIM>& dUdx, const Real rho, const Real delta)
  {
    Real Sijmag = 0.0;
    for (int m = 0; m < AMREX_SPACEDIM; m++) {
      for (int n = 0; n < AMREX_SPACEDIM; n++) {
        Sijmag += (0.5 * (dUdx[m][n] + dUdx[n][m])) *
                  (0.5 * (dUdx[m][n] + dUdx[n][m])); // Sij*Sij
      }
    }
    Sijmag = std::sqrt(2.0 * Sijmag);    
    return (rho * Cs * Cs * delta * delta * Sijmag);
  }

  ///////////////////////

/**
 * @brief Calculate turbulent viscosity using the WALE model.
 *        Note that velocity derivatives are evaluated at cell centre (cc).
 *
 * @param i,j,k     x,y,z index
 * @param q         primitive variables array
 * @param dxinv     cell size, used for calculating velocity derivatives
 * @param delta  filter width
 * @param Cs        Smagorinsky constant, WALE model constant is (10.6)^0.5*Cs
 * @param mu_T[out] output turbulent viscosity divided by Cs^2
 */

/// WALE
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
WALE::mu_T_cc(const int i, const int j, const int k, const Array4<const Real>& q,
              const GpuArray<Real, AMREX_SPACEDIM>& dxinv, const Real delta,
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
  mu_T = q(i, j, k, QRHO) * Cw * Cw * delta * delta * std::pow(DijDij, 1.5) /
         (std::pow(SijSij, 2.5) + std::pow(DijDij, 1.25) +
          std::numeric_limits<Real>::denorm_min());
}



                                                                                          111,1         Bot


///////////////////////////////////////////////////

};

////////////////////////////////////////////////////////////////////////////////


