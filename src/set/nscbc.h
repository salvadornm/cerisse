// NSBC
#include <AMReX_Array4.H>
#include <AMReX_REAL.H>


// Parameters for NSBC conditions are read from input
// relaxation paremeters, type, etc.
// the NSBC/farfiekd conditions follow the process
// 1) compute L from last inner cells
// 2) modify  L to adjust to desired BC
// 3) recast d as a function of L
// 4) add rhs as a function of d

//////////////////////////////////////////////////////////////////////////////////////////////
// compute one-side derivatives (2nd order) 
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void one_side_deriv2(
  amrex::IntVect iv, int idir, int isign,
  amrex::Real deltainv, amrex::Real& dp, amrex::Real& dun, amrex::Real& dut,
  amrex::Real& dutt, amrex::Real& drho, const amrex::Array4<const amrex::Real>& q)
{
  const auto iv_dir = amrex::IntVect::TheDimensionVector(idir) * isign;

  const amrex::Real coef= amrex::Real(isign)*deltainv;

  AMREX_D_TERM(int QUN = QU + idir;, int QUT = QU + (idir + 1) % (AMREX_SPACEDIM - 1);
              ,int QUTT = QU + (idir + 2) % (AMREX_SPACEDIM - 1);)

  // pressure            
  dp = coef * (-1.5 * q(iv, QPRES) + 2.0 * q(iv + iv_dir, QPRES) - 0.5 * q(iv + 2 * iv_dir, QPRES)) ;
  // velocity normal and tangential
  AMREX_D_TERM(dun  = coef * (-1.5 * q(iv, QUN)  + 2.0 * q(iv + iv_dir, QUN)  - 0.5 * q(iv + 2 * iv_dir, QUN));
             , dut  = coef * (-1.5 * q(iv, QUT)  + 2.0 * q(iv + iv_dir, QUT)  - 0.5 * q(iv + 2 * iv_dir, QUT));
             , dutt = coef * (-1.5 * q(iv, QUTT) + 2.0 * q(iv + iv_dir, QUTT) - 0.5 * q(iv + 2 * iv_dir, QUTT)); )
  // density
  drho =  coef* (-1.5 * q(iv, QRHO) + 2.0 * q(iv + iv_dir, QRHO) - 0.5 * q(iv + 2 * iv_dir, QRHO));

#if NUM_SPECIES > 1
  // species
  for (int ns = 0; ns < NUM_SPECIES; ++ns) {
    dY[ns] =  coef*(-1.5 * q(iv, QFS+ns) + 2.0 * q(iv + iv_dir, QFS+ns) - 0.5 * q(iv + 2 * iv_dir, QFS+ns));
  }
#endif

}
//////////////////////////////////////////////////////////////////////////////////////////////
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void one_side_deriv1(
  amrex::IntVect iv, int idir, int isign,
  amrex::Real deltainv, amrex::Real& dp, amrex::Real& dun, amrex::Real& dut,
  amrex::Real& dutt, amrex::Real& drho, amrex::Real* dY, const amrex::Array4<const amrex::Real>& q)
{
  const auto iv_dir = amrex::IntVect::TheDimensionVector(idir) * isign;

  const amrex::Real coef= amrex::Real(isign)*deltainv;

  AMREX_D_TERM(int QUN = QU + idir;, int QUT = QU + (idir + 1) % (AMREX_SPACEDIM - 1);
              ,int QUTT = QU + (idir + 2) % (AMREX_SPACEDIM - 1);)

  // pressure            
  dp = coef * (- q(iv, QPRES) +   q(iv + iv_dir, QPRES)) ;
  // velocity normal and tangential
  AMREX_D_TERM(dun  = coef * (- q(iv, QUN)  +  q(iv + iv_dir, QUN) );
             , dut  = coef * (- q(iv, QUT)  +  q(iv + iv_dir, QUT) );
             , dutt = coef * (- q(iv, QUTT) +  q(iv + iv_dir, QUTT)); )
  // density
  drho =  coef* (- q(iv, QRHO) +  q(iv + iv_dir, QRHO));

#if NUM_SPECIES > 1
  // species
  for (int ns = 0; ns < NUM_SPECIES; ++ns) {
    dY[ns] =  coef* (- q(iv, QFS+ns) +  q(iv + iv_dir, QFS+ns));
  }
#endif

}
//////////////////////////////////////////////////////////////////////////////////////////////
// @brief compute L-waves from q and derivatives (in x)
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void compute_Lwaves_x(
  amrex::IntVect iv, int idir, int isign,
  amrex::Real& dp, amrex::Real& dun, amrex::Real& dut,
  amrex::Real& dutt, amrex::Real& drho, amrex::Real* /*dY*/, const amrex::Array4<const amrex::Real>& q, amrex::Real& L)
  {

    int QUN = QU + idir;
    // velocity and sound spped
    const amrex::Real u   = q(iv,QUN);
    const amrex::Real c   = q(iv,QC);
    const amrex::Real rho = q(iv,QRHO);
    
    L[0] = (u - c)*(dp - rho*c*dun);
    
    // TODO generalise
    L[1] = u*(c*c*drho - dp);   // normal          
    L[2] = u*dut;               // tangential  (y)  2D /3D
    L[3] = u*dutt;              // tangential2 (z)  3D

    L[4] = (u + c)*(dp  + rho*c*dun);
#if NUM_SPECIES > 1
    for (int ns = 0; ns < NUM_SPECIES; ++ns) {
      L[5 + ns] =  u*dY[ns];
    }
#endif

  }
//////////////////////////////////////////////////////////////////////////////////////////////
// @brief compute T-waves from q and tangential derivatives
// TODO 
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void compute_Twaves_x(
  amrex::IntVect iv, int idir, int isign,
  amrex::Real deltainv, amrex::Real& dp, amrex::Real& dun, amrex::Real& dut,
  amrex::Real& dutt, amrex::Real& drho, amrex::Real* dY, const amrex::Array4<const amrex::Real>& q, amrex::Real& T)
  {

    int QUN = QU + idir;
    // velocity and sound spped
    const amrex::Real u   = q(iv,QUN);
    const amrex::Real c   = q(iv,QC);
    const amrex::Real rho = q(iv,QRHO);
    
    T[0] = (u - c)*(dp - rho*c*dun);
    
    // TODO generalise
    T[1] = u*(c*c*drho - dp);
    T[2] = u*dut;
    T[3] = u*dutt;

    T[4] = (u + c)*(dp  + rho*c*dun);
#if NUM_SPECIES > 1
    for (int ns = 0; ns < NUM_SPECIES; ++ns) {
      T[5 + ns] =  u*dY[ns];
    }
#endif

  }
//////////////////////////////////////////////////////////////////////////////////////////////
// @brief compute d-waves L
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void compute_dfromL(
  amrex::IntVect iv, int idir, int isign, const amrex::Array4<const amrex::Real>& q,
  amrex::Real* L, amrex::Real* d)
  {

    const Real rho   = q(iv,QRHO);
    const Real o_rho = 1.0/q(iv,QRHO);

    d[0] =  -L[0] + rho*L[2]/Temp - L[4];   // drho/dt
    d[1] =  -c*(L[5] - L[1])*o_rho;         // dun/dt
    d[2] =  -L[2];                         // dut/dt
    d[3] =  -L[3];                         // dutt/dt
    d[4] =  -coef*(L[0] - L[5])*o_rho - L[1];  // dT/dt
#if NUM_SPECIES > 1
    for (int ns = 0; ns < NUM_SPECIES; ++ns) {
     d[5 + ns] = - L[5+ns];
     d[0] += rho*ratio*L[5+ns];
    }
#endif  

  }

  
  // compute 
  /**
  * @brief Compute rhs in the ghost points based on outflow with target pressure
  * @param Ptarget   target pressure in the outflow
  * @param prims     array of conservative variables
  * @param rhs       rhs
  **/
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE void nscbc_outflow(amrex::IntVect iv,int idir, int isign,
    const cls_t* cls, const Real Ptarget, const Array4<Real>& prims, const Array4<Real>& rhs)
    {

    // calc derivatives
    Real dun, dut,dutt,drho,dp;
    Real dY[NUM_SPECIES];
    one_side_deriv1(iv, idir, isign, deltainv, dp, dun, dut,dutt,rho, dY, prims);
    // calc Lwaves
    Real L[cls.NUM_WAVES];
    compute_Lwaves_x(iv, idir, isign, dp,  dun, dut,dutt, drho, dY, prims, L);

    // modify Lwave
    const int IN_ACOUSTIC = isign; //TODO general based 

    L[IN_ACOUSTIC] = 0.;

    Real d[cls.NUM_WAVES];
    compute_dfromL(iv, idir, isign,prims,L,d);


    // adjust RHS
    rhs(iv,URHO) += d[0];
    rhs(iv,UMX)  += prims(iv,cls.QU)*d[0] + prims(iv,cls.QRHO)*d[1];
    rhs(iv,UMY)  += prims(iv,cls.QV)*d[0] + prims(iv,cls.QRHO)*d[2];
    rhs(iv,UMZ)  += prims(iv,cls.QW)*d[0] + prims(iv,cls.QRHO)*d[3];

  

    // ..



    }



