#ifndef LES_H_
#define LES_H_

#include <CNSconstants.h>

#include "diff_ops.H"

//////////////////////////////// LES TEMPLATE /////////////////////////////////
template <typename param, typename idx_t>
class LES_t {

  public:
  AMREX_GPU_HOST_DEVICE
  LES_t() {}

  AMREX_GPU_HOST_DEVICE
  ~LES_t() {}

  static constexpr Real Pr_o_Prsgs = param::Pr_o_Prsgs;
  const int order  = param::order;
  static constexpr Real Scsgs_inv = 1.0/param::Scsgs;
  // WALE constant
  //static constexpr std::sqrt(10.6)=> 3.255764119219941 
  static constexpr Real Cw = 3.255764119219941 * param::Cs;
    
  /**
   * \brief calculates filter width
   * \param[in] dx
   * \return Delta
  */
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE  
  Real calc_delta(  const GpuArray<Real, AMREX_SPACEDIM>& dx) const{

    if constexpr (param::fixDelta) {
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
   * \brief calculates sub-grid bulk viscosity at FACE i+1/2??? CHECK 
   *  cell centred
   * \param[in]  i,j,k
   * \param[in]  q     primitive variables array 
   * \param  dxinv     cell size, used for calculating velocity derivatives
   * \param[in]  Delta
   * \param[out] xi_sgs
  */
  // central in face (depending on idir)
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE Real xi_sgs_face(
  const int i, const int j, const int k, const int d1,  const Array4<const Real>& q,
  const GpuArray<Real, AMREX_SPACEDIM>& dxinv, const Real delta)
  {
    // Calculate derivatives at cell centers using  central differences at cell faces
    const amrex::IntVect iv{AMREX_D_DECL(i, j, k)};
    Real dUdx[3][3] = {{0.0}};

    // directions
    const int d2 = d1 == 0 ? 1 : 0;
    const int d3 = d1 == 2 ? 1 : 2;
    AMREX_D_TERM(const int QU1 = idx_t::QU + d1;, const int QU2 = idx_t::QU + d2;
               , const int QU3 = idx_t::QU + d3;) 

    dUdx[0][0]  = normal_diff<param::order>(iv, d1, QU1, q, dxinv);   // dudx00
#if (AMREX_SPACEDIM >= 2)
    dUdx[1][0]  = normal_diff<param::order>(iv, d1, QU2, q, dxinv);
    dUdx[0][1]  = tangent_diff<param::order>(iv, d1, d2, QU1, q, dxinv);
    dUdx[1][1]  = tangent_diff<param::order>(iv, d1, d2, QU2, q, dxinv);
#endif
#if (AMREX_SPACEDIM == 3)
    dUdx[2][0]  = normal_diff<param::order>(iv, d1, QU3, q, dxinv);
    dUdx[0][2]  = tangent_diff<param::order>(iv, d1, d3, QU1, q, dxinv);
    dUdx[2][2]  = tangent_diff<param::order>(iv, d1, d3, QU3, q, dxinv);
#endif  

    // || Sij ||
    Real Sijmag = 0.0;
    for (int m = 0; m < AMREX_SPACEDIM; m++) {
      for (int n = 0; n < AMREX_SPACEDIM; n++) {
        Sijmag += (0.5 * (dUdx[m][n] + dUdx[n][m])) *
                  (0.5 * (dUdx[m][n] + dUdx[n][m])); // Sij*Sij
      }
    }
    //Sij2 = ||Sij||**2
    const Real Sij2 = 2.0*Sijmag;
    return( two_third*q(i, j, k, idx_t::QRHO) *param::CI  * delta * delta *Sij2);
  }
  /**
   * \brief calculates sub-grid velocity at cell using vorticity (2nd order only)
   *  
   * \param[in]  i,j,k
   * \param[in]  q     primitive variables array 
   * \param  dxinv     cell size, used for calculating velocity derivatives
   * \param[in]  Delta
   * \param[out] usgs   sub-grid velocity
  */
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE Real usgs_cell(
  const int i, const int j, const int k, const Array4<const Real>& q,
  const GpuArray<Real, AMREX_SPACEDIM>& dxinv, const Real Delta) const
  {
    // Unit vectors in index space (x,y,z)
    constexpr int ex = 0, ey = 1, ez = 2;

    // Helper: compute vorticity vector at a given cell iv from velocity gradients
    auto vort_at = [&](const IntVect& iv) -> amrex::GpuArray<Real,3>
    {
      Real dUdx[3][3] = {{0.0_rt}}; // dUdx[m][n] = d(u_m)/d(x_n)

      for (int m = 0; m < AMREX_SPACEDIM; ++m) {
        for (int n = 0; n < AMREX_SPACEDIM; ++n) {
          dUdx[m][n] = normal_diff_cc<param::order>(iv, n, idx_t::QU + m, q, dxinv);
        }
      }

      amrex::GpuArray<Real,3> w{0.0_rt, 0.0_rt, 0.0_rt};

#if (AMREX_SPACEDIM == 3)
      // ω = ∇×u
      // ωx = ∂w/∂y - ∂v/∂z
      w[0] = dUdx[2][ey] - dUdx[1][ez];
      // ωy = ∂u/∂z - ∂w/∂x
      w[1] = dUdx[0][ez] - dUdx[2][ex];
      // ωz = ∂v/∂x - ∂u/∂y
      w[2] = dUdx[1][ex] - dUdx[0][ey];
#else
      // 2D: only ωz is non-zero (out-of-plane vorticity)
      w[2] = dUdx[1][ex] - dUdx[0][ey];
#endif
      return w;
    };

    const IntVect iv{AMREX_D_DECL(i,j,k)};

    // Vorticity at center
    const auto w0 = vort_at(iv);

    // Laplacian of vorticity: ∇²ω ≈ Σ_n (ω(i+1)-2ω(i)+ω(i-1)) * (1/Δx_n^2)
    amrex::GpuArray<Real,3> lapw{0.0_rt, 0.0_rt, 0.0_rt};

    // x-direction contribution
    {
      const auto wp = vort_at(iv + IntVect{AMREX_D_DECL(1,0,0)});
      const auto wm = vort_at(iv + IntVect{AMREX_D_DECL(-1,0,0)});
      const Real c = dxinv[0]*dxinv[0];
      lapw[0] += (wp[0] - 2.0_rt*w0[0] + wm[0]) * c;
      lapw[1] += (wp[1] - 2.0_rt*w0[1] + wm[1]) * c;
      lapw[2] += (wp[2] - 2.0_rt*w0[2] + wm[2]) * c;
    }
    
    // y-direction contribution
    {
      const auto wp = vort_at(iv + IntVect{AMREX_D_DECL(0,1,0)});
      const auto wm = vort_at(iv + IntVect{AMREX_D_DECL(0,-1,0)});
      const Real c = dxinv[1]*dxinv[1];
      lapw[0] += (wp[0] - 2.0_rt*w0[0] + wm[0]) * c;
      lapw[1] += (wp[1] - 2.0_rt*w0[1] + wm[1]) * c;
      lapw[2] += (wp[2] - 2.0_rt*w0[2] + wm[2]) * c;
    }

#if (AMREX_SPACEDIM == 3)
  // z-direction contribution
    {
      const auto wp = vort_at(iv + IntVect{AMREX_D_DECL(0,0,1)});
      const auto wm = vort_at(iv + IntVect{AMREX_D_DECL(0,0,-1)});
      const Real c = dxinv[2]*dxinv[2];
      lapw[0] += (wp[0] - 2.0_rt*w0[0] + wm[0]) * c;
      lapw[1] += (wp[1] - 2.0_rt*w0[1] + wm[1]) * c;
      lapw[2] += (wp[2] - 2.0_rt*w0[2] + wm[2]) * c;
    }
#endif

    // |∇²(∇×u)|
#if (AMREX_SPACEDIM == 3)
    const Real mag = std::sqrt(lapw[0]*lapw[0] + lapw[1]*lapw[1] + lapw[2]*lapw[2]);
#else
    const Real mag = std::abs(lapw[2]); // only z-component exists in 2D
#endif

    constexpr Real C2 = 0.1;

    // u' = C2 Δ^3 * | ... |

    return (C2 * (Delta*Delta*Delta) * mag);
  }


////////////////////////////////////////////////////////////////
};

//////////////////////////////// SMAG TEMPLATE /////////////////////////////////
template <typename param, typename idx_t>
class Smagorinsky_t : public LES_t<param, idx_t> {

  public :

  /**
   * \brief calculates sub-grid viscosity
   * \param[in] i,j,k,dxinv,delta
   * \param[out] mu_T 
  */
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE void visc_sgs(
  const int i, const int j, const int k, const Array4<const Real>& q,
  const GpuArray<Real, AMREX_SPACEDIM>& dxinv, const Real delta, Real& mu_T) const
  {
    // Calculate derivatives at cell centers uisng  central differences
    const amrex::IntVect iv{AMREX_D_DECL(i, j, k)};
    Real dUdx[3][3] = {{0.0}};
    // finite difference central order 2/4/6
    for (int m = 0; m < AMREX_SPACEDIM; m++) {
      for (int n = 0; n < AMREX_SPACEDIM; n++) {
      	dUdx[m][n] = normal_diff_cc<param::order>(iv, n, idx_t::QU + m, q, dxinv);
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
    mu_T = q(i, j, k, idx_t::QRHO) * param::Cs * param::Cs * delta * delta * Sijmag;
  }
  /**
   * \brief calculates sub-grid conductivity
   * \param[in] i,j,k,dxinv,delta
   * \param[out] cond_T 
  */
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE void cond_sgs(
  const int i, const int j, const int k, const Array4<const Real>& q,
  const GpuArray<Real, AMREX_SPACEDIM>& dxinv, const Real delta, const Real& Cp_o_Pr, Real& cond_T)  const
  {
    Real mu_T;
    visc_sgs(i,j,k,q,dxinv,delta,mu_T);
    cond_T = mu_T*this->Pr_o_Prsgs*Cp_o_Pr;
  }  
  /**
   * \brief calculates sub-grid diffusivity
   * \param[in] i,j,k,dxinv,delta
   * \param[out]  rhoD_T 
  */
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE void diff_sgs(
  const int i, const int j, const int k, const Array4<const Real>& q,
  const GpuArray<Real, AMREX_SPACEDIM>& dxinv,const Real delta, Real& rhoD_T) const
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
  const GpuArray<Real, AMREX_SPACEDIM>& dxinv, const Real delta, const Real& Cp_o_Pr, 
  Real& mu_T, Real& cond_T, Real& rhoD_T) const
  {
    visc_sgs(i,j,k,q,dxinv,delta,mu_T);
    cond_T = mu_T*param::Pr_o_Prsgs*Cp_o_Pr;
    rhoD_T = mu_T*this->Scsgs_inv;
  }  
  /**
   * \brief calculates sub-grid time scale
   * \param[in] i,j,k,dxinv,delta
   * \param[out]  tau_T 
  */
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE Real tau_sgs(
  const int i, const int j, const int k, const Array4<const Real>& q,
  const GpuArray<Real, AMREX_SPACEDIM>& dxinv,const Real delta) const
  {
    Real mu_T;
    visc_sgs(i,j,k,q,dxinv,delta,mu_T);
    Real tau_T = q(i, j, k, idx_t::QRHO)*delta*delta/(mu_T + smalleps);
    return tau_T;
  }  

};
//////////////////////////////// WALE TEMPLATE /////////////////////////////////
template <typename param, typename idx_t>
class WALE_t : public LES_t<param, idx_t> {

  using LES_t<param,idx_t>::Cw;

public :

  AMREX_GPU_DEVICE AMREX_FORCE_INLINE void visc_sgs(
  const int i, const int j, const int k, const Array4<const Real>& q,
  const GpuArray<Real, AMREX_SPACEDIM>& dxinv, const Real delta, Real& mu_T) const
  {
    // Calculate derivatives at cell centers, second order central difference
    const amrex::IntVect iv{AMREX_D_DECL(i, j, k)};
    Real dUdx[3][3] = {{0.0}};
    // finite difference central order 2/4/6
    for (int m = 0; m < AMREX_SPACEDIM; m++) {
      for (int n = 0; n < AMREX_SPACEDIM; n++) {
        //dUdx[m][n]  = normal_diff_cc<param::order>(iv, n, this->QUn[m], q, dxinv); // dUmdn
	dUdx[m][n] = normal_diff_cc<param::order>(iv, n, idx_t::QU + m, q, dxinv);

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

    mu_T = q(i, j, k, idx_t::QRHO) * Cw * Cw * delta * delta * std::pow(DijDij, 1.5) /
         (std::pow(SijSij, 2.5) + std::pow(DijDij, 1.25) +
          std::numeric_limits<Real>::denorm_min());
  }
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE void cond_sgs(
  const int i, const int j, const int k, const Array4<const Real>& q,
  const GpuArray<Real, AMREX_SPACEDIM>& dxinv, const Real delta, const Real& Cp_o_Pr, Real& cond_T)
  {
    Real mu_T;
    visc_sgs(i,j,k,q,dxinv,delta,mu_T);
    cond_T = mu_T*this->Pr_o_Prsgs*Cp_o_Pr;
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
  const GpuArray<Real, AMREX_SPACEDIM>& dxinv, const Real delta, const Real& Cp_o_Pr, 
  Real& mu_T, Real& cond_T, Real& rhoD_T) const
  {
    visc_sgs(i,j,k,q,dxinv,delta,mu_T);
    cond_T = mu_T*this->Pr_o_Prsgs*Cp_o_Pr;
    rhoD_T = mu_T*this->Scsgs_inv;
  }  
  /**
   * \brief calculates sub-grid time scale
   * \param[in] i,j,k,dxinv,delta
   * \param[out]  tau_T 
  */
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE Real tau_sgs(
  const int i, const int j, const int k, const Array4<const Real>& q,
  const GpuArray<Real, AMREX_SPACEDIM>& dxinv,const Real delta) const
  {
    Real mu_T;
    visc_sgs(i,j,k,q,dxinv,delta,mu_T);
    Real tau_T = q(i, j, k, this->QRHO)*delta*delta/(mu_T + smalleps);
    return tau_T;
  }  
};

////////////////////////////////////////////////////////////////////////////////
#endif

