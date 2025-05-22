#ifndef LES_H_
#define LES_H_

#include <CNSconstants.h>

//////////////////////////////// LES TEMPLATE /////////////////////////////////
template <typename param, typename cls_t>
class LES_t {

  public:
  AMREX_GPU_HOST_DEVICE
  LES_t() {}

  AMREX_GPU_HOST_DEVICE
  ~LES_t() {}

  const Real Prsgs = param::Prsgs;
  const Real Scsgs = param::Scsgs;
  const Real Delta = param::Delta;
  const bool fixDelta = param::fixDelta;
  const Real CI = param::CI;  // Yoshizawa constant
  const Real Cs = param::Cs;  // Smagorinky constant

  const int order  = param::order;
  static constexpr Real Prsgs_inv = 1.0/param::Prsgs;
  static constexpr Real Scsgs_inv = 1.0/param::Scsgs;
  // Indexes
  static constexpr int QUn[3]={cls_t::QU,cls_t::QV,cls_t::QW};
  // WALE constant
  static constexpr Real Cw = std::sqrt(10.6) * param::Cs;
    
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
template <typename param, typename cls_t>
class Smagorinsky_t : public LES_t<param, cls_t> {

  public :

  /**
   * \brief calculates sub-grid viscosiy
   * \param[in] i,j,k,dxinv,delta
   * \param[out] mu_T 
  */
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE void visc_sgs(
  const int i, const int j, const int k, const Array4<const Real>& q,
  const GpuArray<Real, AMREX_SPACEDIM>& dxinv, const Real delta, Real& mu_T)
  {
    // Calculate derivatives at cell centers, second order central difference
    const amrex::IntVect iv{AMREX_D_DECL(i, j, k)};
    Real dUdx[3][3] = {{0.0}};
    // finite difference central order 2/4/6
    for (int m = 0; m < AMREX_SPACEDIM; m++) {
      for (int n = 0; n < AMREX_SPACEDIM; n++) {
        dUdx[m][n]  = normal_diff<param::order>(iv, n, this->QUn[m], q, dxinv); // dUmdn
     }
    }
    // || Sij ||
    Real Sijmag = 0.0;
    for (int m = 0; m < AMREX_SPACEDIM; m++) {
      for (int n = 0; n < AMREX_SPACEDIM; n++) {
        Sijmag += (0.5 * (dUdx[m][n] + dUdx[n][m])) *
                (0.5 * (dUdx[m][n] + dUdx[n][m])); // Sij*Sij
      }
    }
    Sijmag = std::sqrt(2.0 * Sijmag);
    mu_T = q(i, j, k, this->QRHO) * param::Cs * param::Cs * delta * delta * Sijmag;
  }
  /**
   * \brief calculates sub-grid conductivity
   * \param[in] i,j,k,dxinv,delta
   * \param[out] cond_T 
  */
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE void cond_sgs(
  const int i, const int j, const int k, const Array4<const Real>& q,
  const GpuArray<Real, AMREX_SPACEDIM>& dxinv, const Real delta, const Real& Cp, Real& cond_T)
  {
    Real mu_T;
    visc_sgs(i,j,k,q,dxinv,delta,mu_T);
    cond_T = mu_T*this->Prsgs_inv*Cp;
  }  
  /**
   * \brief calculates sub-grid diffusivity
   * \param[in] i,j,k,dxinv,delta
   * \param[out]  rhoD_T 
  */
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE void diff_sgs(
  const int i, const int j, const int k, const Array4<const Real>& q,
  const GpuArray<Real, AMREX_SPACEDIM>& dxinv,const Real delta, Real& rhoD_T)
  {
    Real mu_T;
    visc_sgs(i,j,k,q,dxinv,delta,mu_T);
    rhoD_T = mu_T*this->Scsgs_inv;
  }  
  /**
   * \brief calculates all sgs properties, viscosity, conductivity and diffusivity
   * \param[in] i,j,k,dxinv,delta
   * \param[out]  rhoD_T 
  */
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE void compute_sgsterms(
  const int i, const int j, const int k, const Array4<const Real>& q,
  const GpuArray<Real, AMREX_SPACEDIM>& dxinv, const Real delta, const Real& Cp, 
  Real& mu_T, Real& cond_T, Real& rhoD_T)
  {
    visc_sgs(i,j,k,q,dxinv,delta,mu_T);
    cond_T = mu_T*this->Prsgs_inv*Cp;
    rhoD_T = mu_T*this->Scsgs_inv;
  }  
  /**
   * \brief calculates sub-grid time scale
   * \param[in] i,j,k,dxinv,delta
   * \param[out]  tau_T 
  */
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE void tau_sgs(
  const int i, const int j, const int k, const Array4<const Real>& q,
  const GpuArray<Real, AMREX_SPACEDIM>& dxinv,const Real delta, Real& tau_T)
  {
    Real mu_T;
    visc_sgs(i,j,k,q,dxinv,delta,mu_T);
    tau_T = q(i, j, k, this->QRHO)*delta*delta/(mu_T + smalleps);
  }  

};
//////////////////////////////// SMAG TEMPLATE /////////////////////////////////
template <typename param, typename cls_t>
class WALE_t : public LES_t<param, cls_t> {

public :

  AMREX_GPU_DEVICE AMREX_FORCE_INLINE void visc_sgs(
  const int i, const int j, const int k, const Array4<const Real>& q,
  const GpuArray<Real, AMREX_SPACEDIM>& dxinv, const Real delta, Real& mu_T)
  {
    // Calculate derivatives at cell centers, second order central difference
    const amrex::IntVect iv{AMREX_D_DECL(i, j, k)};
    Real dUdx[3][3] = {{0.0}};
    // finite difference central order 2/4/6
    for (int m = 0; m < AMREX_SPACEDIM; m++) {
      for (int n = 0; n < AMREX_SPACEDIM; n++) {
        dUdx[m][n]  = normal_diff<param::order>(iv, n, this->QUn[m], q, dxinv); // dUmdn
     }
    }
    //
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
    Real DijDij = 0.0; Real Dkk = divu * divu * one_third; 
    for (int m = 0; m < AMREX_SPACEDIM; m++) {
      for (int n = 0; n < AMREX_SPACEDIM; n++) {
        SijSij += (0.5 * (dUdx[m][n] + dUdx[n][m])) *
                (0.5 * (dUdx[m][n] + dUdx[n][m])); // Sij*Sij
        DijDij +=
          (0.5 * (dUdx2[m][n] + dUdx2[n][m]) - Real(m == n) * Dkk) *
          (0.5 * (dUdx2[m][n] + dUdx2[n][m]) - Real(m == n) * Dkk); // Dij*Dij
      }
    }

    mu_T = q(i, j, k, this->QRHO) * param::Cw * param::Cw * delta * delta * std::pow(DijDij, 1.5) /
         (std::pow(SijSij, 2.5) + std::pow(DijDij, 1.25) +
          std::numeric_limits<Real>::denorm_min());
  }
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE void cond_sgs(
  const int i, const int j, const int k, const Array4<const Real>& q,
  const GpuArray<Real, AMREX_SPACEDIM>& dxinv, const Real delta, const Real& Cp, Real& cond_T)
  {
    Real mu_T;
    visc_sgs(i,j,k,q,dxinv,delta,mu_T);
    cond_T = mu_T*this->Prsgs_inv*Cp;
  }  
  /**
   * \brief calculates sub-grid diffusivity
   * \param[in] i,j,k,dxinv,delta
   * \param[out]  rhoD_T 
  */
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE void diff_sgs(
  const int i, const int j, const int k, const Array4<const Real>& q,
  const GpuArray<Real, AMREX_SPACEDIM>& dxinv,const Real delta, Real& rhoD_T)
  {
    Real mu_T;
    visc_sgs(i,j,k,q,dxinv,delta,mu_T);
    rhoD_T = mu_T*this->Scsgs_inv;
  }  
  /**
   * \brief calculates all sgs properties, viscosity, conductivity and diffusivity
   * \param[in] i,j,k,dxinv,delta
   * \param[out]  rhoD_T 
  */
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE void compute_sgsterms(
  const int i, const int j, const int k, const Array4<const Real>& q,
  const GpuArray<Real, AMREX_SPACEDIM>& dxinv, const Real delta, const Real& Cp, 
  Real& mu_T, Real& cond_T, Real& rhoD_T)
  {
    visc_sgs(i,j,k,q,dxinv,delta,mu_T);
    cond_T = mu_T*this->Prsgs_inv*Cp;
    rhoD_T = mu_T*this->Scsgs_inv;
  }  
  /**
   * \brief calculates sub-grid time scale
   * \param[in] i,j,k,dxinv,delta
   * \param[out]  tau_T 
  */
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE void tau_sgs(
  const int i, const int j, const int k, const Array4<const Real>& q,
  const GpuArray<Real, AMREX_SPACEDIM>& dxinv,const Real delta, Real& tau_T)
  {
    Real mu_T;
    visc_sgs(i,j,k,q,dxinv,delta,mu_T);
    tau_T = q(i, j, k, this->QRHO)*delta*delta/(mu_T + smalleps);
  }  
};

////////////////////////////////////////////////////////////////////////////////
#endif

