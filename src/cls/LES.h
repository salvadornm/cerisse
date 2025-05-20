#ifndef LES_H_
#define LES_H_

//////////////////////////////// LES TEMPLATE /////////////////////////////////
template <typename param, typename cls_t>
class LES_t {
  private:
  Real Prsgs = param::Prsgs;
  Real Scsgs = param::Scsgs;
  Real Delta = param::Delta;
  bool fixDelta = param::fixDelta;
  Real CI = param::CI;  // Yoshizawa constant


  static const int QU=cls_t::QU;
  static const int QV=cls_t::QV;
  static const int QW=cls_t::QW;
  
  public:
  
  /**
   * \brief calculates filter width
   * \param[in] dx
   * \return Delta
  */
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE Real calc_delta(
  const GpuArray<Real, AMREX_SPACEDIM>& dx){

    if (param::fixDelta) {
      return(param::Delta);
    }
    else
    {
#if (AMREX_SPACEDIM==1)
      return(dx[0]);
#elif (AMREX_SPACEDIM==2)
      return(std::sqrt(dx[0]*dx[1]));
#else    
      return(std::cbrt(dx[0]*dx[1]*dx[2]));
#endif    
    }
  }

  /**
   * \brief calculates velocity derivatives and store in dUdx[n][m]
   *  cell centred
   * \param[in] dx
   * \return Delta
  */
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE void gradvel(
  const int i, const int j, const int k, const Array4<const Real>& q,
  const GpuArray<Real, AMREX_SPACEDIM>& dxinv, const GpuArray<Real, AMREX_SPACEDIM>& dUdx)
  {

    dUdx[0][0] = 0.5 * dxinv[0] * (q(i + 1, j, k, QU) - q(i - 1, j, k, QU)); // dudx
# if AMREX_SPACEDIM > 1    
    dUdx[1][0] = 0.5 * dxinv[0] * (q(i + 1, j, k, QV) - q(i - 1, j, k, QV)); // dvdx
    dUdx[0][1] = 0.5 * dxinv[1] * (q(i, j + 1, k, QU) - q(i, j - 1, k, QU)); // dudy
    dUdx[1][1] = 0.5 * dxinv[1] * (q(i, j + 1, k, QV) - q(i, j - 1, k, QV)); // dvdy
#endif    
#if AMREX_SPACEDIM == 3
    dUdx[2][0] = 0.5 * dxinv[0] * (q(i + 1, j, k, QW) - q(i - 1, j, k, QW)); // dwdx
    dUdx[2][1] = 0.5 * dxinv[1] * (q(i, j + 1, k, QW) - q(i, j - 1, k, QW)); // dwdy
    dUdx[0][2] = 0.5 * dxinv[2] * (q(i, j, k + 1, QU) - q(i, j, k - 1, QU)); // dudz
    dUdx[1][2] = 0.5 * dxinv[2] * (q(i, j, k + 1, QV) - q(i, j, k - 1, QV)); // dvdz
    dUdx[2][2] = 0.5 * dxinv[2] * (q(i, j, k + 1, QW) - q(i, j, k - 1, QW)); // dwdz
#endif
  }

  /**
   * \brief calculates sub-grid bulk viscoisty 
   *  cell centred
   * \param[in]  i,j,k
   * \param[in]  q     primitive variables array 
   * \param  dxinv     cell size, used for calculating velocity derivatives
   * \param[in]  Delta
   * \param[out] xi_sgs
  */



////////////////////////////////////////////////////////////////
};

//////////////////////////////// SMAG TEMPLATE /////////////////////////////////
// inhertance from LES


////////////////////////////////////////////////////////////////////////////////
#endif

