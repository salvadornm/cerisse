#ifndef CNS_DIFFUSION_K_H_
#define CNS_DIFFUSION_K_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_CONSTANTS.H>
#include <cmath>

#include "index_macros.H"

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void
cns_diff_x (int i, int j, int k,
            amrex::Array4<amrex::Real const> const& q,
            amrex::Array4<amrex::Real const> const& coeffs,
            amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dxinv,
            amrex::Array4<amrex::Real> const& fx) noexcept
{
  using amrex::Real;

  Real dTdx = (q(i,j,k,QTEMP)-q(i-1,j,k,QTEMP))*dxinv[0];
  AMREX_D_TERM(Real  dudx = (q(i,j,k,QU)-q(i-1,j,k,QU))*dxinv[0];,
               Real  dvdx = (q(i,j,k,QV)-q(i-1,j,k,QV))*dxinv[0];,
               Real  dwdx = (q(i,j,k,QW)-q(i-1,j,k,QW))*dxinv[0];);
#if (AMREX_SPACEDIM >= 2)
  Real dudy = (q(i,j+1,k,QU)+q(i-1,j+1,k,QU)-q(i,j-1,k,QU)-q(i-1,j-1,k,QU))*(0.25*dxinv[1]);
  Real dvdy = (q(i,j+1,k,QV)+q(i-1,j+1,k,QV)-q(i,j-1,k,QV)-q(i-1,j-1,k,QV))*(0.25*dxinv[1]);
#endif
#if (AMREX_SPACEDIM == 3)
  Real dudz = (q(i,j,k+1,QU)+q(i-1,j,k+1,QU)-q(i,j,k-1,QU)-q(i-1,j,k-1,QU))*(0.25*dxinv[2]);
  Real dwdz = (q(i,j,k+1,QW)+q(i-1,j,k+1,QW)-q(i,j,k-1,QW)-q(i-1,j,k-1,QW))*(0.25*dxinv[2]);
#endif
  Real divu = AMREX_D_TERM(dudx, + dvdy, + dwdz);

  Real mu = 0.5*(coeffs(i,j,k,CMU)+coeffs(i-1,j,k,CMU));
  Real xi = 0.5*(coeffs(i,j,k,CXI)+coeffs(i-1,j,k,CXI));
  Real lamda = 0.5*(coeffs(i,j,k,CLAM)+coeffs(i-1,j,k,CLAM));
  Real rhoD[NUM_SPECIES];
  for (int n = 0; n < NUM_SPECIES; ++n)
    rhoD[n] = 0.5*(coeffs(i,j,k,CRHOD+n)+coeffs(i-1,j,k,CRHOD+n));

  if (CNS::do_les) {
#if AMREX_SPACEDIM == 2
    Real dwdx = 0.0;
    Real dwdy = 0.0;
    Real dwdz = 0.0;
    Real dudz = 0.0;
    Real dvdz = 0.0;
    Real delta = std::pow(dxinv[0]*dxinv[1], -1./2.);
#else
    Real dwdy = (q(i,j+1,k,QW)+q(i-1,j+1,k,QW)-q(i,j-1,k,QW)-q(i-1,j-1,k,QW))*(0.25*dxinv[1]);
    Real dvdz = (q(i,j,k+1,QV)+q(i-1,j,k+1,QV)-q(i,j,k-1,QV)-q(i-1,j,k-1,QV))*(0.25*dxinv[2]);
    Real delta = std::pow(dxinv[0]*dxinv[1]*dxinv[2], -1./3.);
#endif
    Real rho = 0.5*(q(i,j,k,QRHO)+q(i-1,j,k,QRHO));
    Real mu_T;

    CNS::les_model->mu_T_fc(dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz,
                            rho,delta,CNS::Cs,mu_T);

    mu += mu_T;

    Real y[NUM_SPECIES];
    for (int n = 0; n < NUM_SPECIES; ++n)
      y[n] = 0.5*(q(i,j,k,QFS+n)+q(i-1,j,k,QFS+n));
    auto eos = pele::physics::PhysicsType::eos();
    Real cp;
    eos.RTY2Cp(rho, 0.5*(q(i,j,k,QTEMP)+q(i-1,j,k,QTEMP)), y, cp);
    lamda += cp*mu_T/CNS::Pr_T;

    for (int n = 0; n < NUM_SPECIES; ++n)
      rhoD[n] += mu_T/CNS::Sc_T;
  }

  AMREX_D_TERM(Real tauxx = mu*(2.0*dudx-2.0/3.0*divu) + xi*divu; ,
               Real tauxy = mu*(dudy+dvdx); ,
               Real tauxz = mu*(dudz+dwdx););

  AMREX_D_TERM(fx(i,j,k,UMX)   += -tauxx; ,
               fx(i,j,k,UMY)   += -tauxy; ,
               fx(i,j,k,UMZ)   += -tauxz;);
  fx(i,j,k,UEDEN) += -0.5*( AMREX_D_TERM( (q(i,j,k,QU)+q(i-1,j,k,QU))*tauxx,
                                         +(q(i,j,k,QV)+q(i-1,j,k,QV))*tauxy,
                                         +(q(i,j,k,QW)+q(i-1,j,k,QW))*tauxz) )
                     -lamda*dTdx;

  // Species transport
  const amrex::IntVect iv{AMREX_D_DECL(i, j, k)};
  const amrex::IntVect ivm(iv - amrex::IntVect::TheDimensionVector(0));

  // Get massfrac, molefrac, enthalpy
  Real mass1[NUM_SPECIES], mass2[NUM_SPECIES];
  Real mole1[NUM_SPECIES], mole2[NUM_SPECIES];
  Real hi1[NUM_SPECIES], hi2[NUM_SPECIES];
  for (int n = 0; n < NUM_SPECIES; ++n) {
    mass1[n] = q(iv,QFS+n);
    mass2[n] = q(ivm,QFS+n);
  }
  auto eos = pele::physics::PhysicsType::eos();
  eos.Y2X(mass1, mole1);
  eos.Y2X(mass2, mole2);

  // Compute species and enthalpy fluxes for ideal EOS
  // Get species/enthalpy diffusion, compute correction vel
  Real T = q(iv,QTEMP);
  eos.T2Hi(T, hi1);
  T = q(ivm,QTEMP);
  eos.T2Hi(T, hi2);

  Real Vc = 0.0;
  const amrex::Real dpdx = (q(iv,QPRES)-q(ivm,QPRES))*dxinv[0];
  const amrex::Real dlnp = dpdx / (0.5*(q(iv,QPRES)+q(ivm,QPRES)));
  for (int n = 0; n < NUM_SPECIES; ++n) {
    const amrex::Real Xface = 0.5*(mole1[n]+mole2[n]);
    const amrex::Real Yface = 0.5*(mass1[n]+mass2[n]);
    const amrex::Real hface = 0.5*(hi1[n]+hi2[n]);
    const amrex::Real dXdx = (mole1[n]-mole2[n])*dxinv[0];
    const amrex::Real Vd = -rhoD[n]*(dXdx+(Xface-Yface)*dlnp);
    Vc += Vd;
    fx(iv,UFS+n) += Vd;
    fx(iv,UEDEN) += Vd*hface;
  }
  // Add correction velocity to fluxes so sum(Vd) = 0
  for (int ns = 0; ns < NUM_SPECIES; ++ns) {
    const amrex::Real Yface = 0.5*(mass1[ns]+mass2[ns]);
    const amrex::Real hface = 0.5*(hi1[ns]+hi2[ns]);
    fx(iv,UFS+ns) -= Yface*Vc;
    fx(iv,UEDEN) -= Yface*hface*Vc;
  }
}

#if (AMREX_SPACEDIM >= 2)
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void
cns_diff_y (int i, int j, int k, amrex::Array4<amrex::Real const> const& q,
            amrex::Array4<amrex::Real const> const& coeffs,
            amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dxinv,
            amrex::Array4<amrex::Real> const& fy) noexcept
{
  using amrex::Real;

  Real  dTdy = (q(i,j,k,QTEMP)-q(i,j-1,k,QTEMP))*dxinv[1];
  AMREX_D_TERM(Real  dudy = (q(i,j,k,QU)-q(i,j-1,k,QU))*dxinv[1];,
               Real  dvdy = (q(i,j,k,QV)-q(i,j-1,k,QV))*dxinv[1];,
               Real  dwdy = (q(i,j,k,QW)-q(i,j-1,k,QW))*dxinv[1];);
  Real  dudx = (q(i+1,j,k,QU)+q(i+1,j-1,k,QU)-q(i-1,j,k,QU)-q(i-1,j-1,k,QU))*(0.25*dxinv[0]);
  Real  dvdx = (q(i+1,j,k,QV)+q(i+1,j-1,k,QV)-q(i-1,j,k,QV)-q(i-1,j-1,k,QV))*(0.25*dxinv[0]);
#if (AMREX_SPACEDIM == 3)
  Real  dvdz = (q(i,j,k+1,QV)+q(i,j-1,k+1,QV)-q(i,j,k-1,QV)-q(i,j-1,k-1,QV))*(0.25*dxinv[2]);
  Real  dwdz = (q(i,j,k+1,QW)+q(i,j-1,k+1,QW)-q(i,j,k-1,QW)-q(i,j-1,k-1,QW))*(0.25*dxinv[2]);
#endif
  Real  divu = AMREX_D_TERM(dudx, + dvdy, + dwdz);

  Real  mu = 0.5*(coeffs(i,j,k,CMU)+coeffs(i,j-1,k,CMU));
  Real  xi = 0.5*(coeffs(i,j,k,CXI)+coeffs(i,j-1,k,CXI));
  Real  lambda = 0.5*(coeffs(i,j,k,CLAM)+coeffs(i,j-1,k,CLAM));
  AMREX_D_TERM(Real tauyy = mu*(2.0*dvdy-(2.0/3.0)*divu) + xi*divu; ,
               Real tauxy = mu*(dudy+dvdx); ,
               Real tauyz = mu*(dwdy+dvdz););

  AMREX_D_TERM(fy(i,j,k,UMX) += -tauxy; ,
                fy(i,j,k,UMY) += -tauyy; ,
                fy(i,j,k,UMZ) += -tauyz;);
  fy(i,j,k,UEDEN) += -0.5*( AMREX_D_TERM( (q(i,j,k,QU)+q(i,j-1,k,QU))*tauxy,
                                         +(q(i,j,k,QV)+q(i,j-1,k,QV))*tauyy,
                                         +(q(i,j,k,QW)+q(i,j-1,k,QW))*tauyz) )
                     -lambda*dTdy;

  // Species transport
  const amrex::IntVect iv{AMREX_D_DECL(i, j, k)};
  const amrex::IntVect ivm(iv - amrex::IntVect::TheDimensionVector(1));

  // Get massfrac, molefrac, enthalpy
  Real mass1[NUM_SPECIES], mass2[NUM_SPECIES];
  Real mole1[NUM_SPECIES], mole2[NUM_SPECIES];
  Real hi1[NUM_SPECIES], hi2[NUM_SPECIES];
  for (int n = 0; n < NUM_SPECIES; ++n) {
    mass1[n] = q(iv,QFS+n);
    mass2[n] = q(ivm,QFS+n);
  }
  auto eos = pele::physics::PhysicsType::eos();
  eos.Y2X(mass1, mole1);
  eos.Y2X(mass2, mole2);

  // Compute species and enthalpy fluxes for ideal EOS
  // Get species/enthalpy diffusion, compute correction vel
  Real T = q(iv,QTEMP);
  eos.T2Hi(T, hi1);
  T = q(ivm,QTEMP);
  eos.T2Hi(T, hi2);

  Real Vc = 0.0;
  const amrex::Real dpdy = (q(iv,QPRES)-q(ivm,QPRES))*dxinv[1];
  const amrex::Real dlnp = dpdy / (0.5*(q(iv,QPRES)+q(ivm,QPRES)));
  for (int n = 0; n < NUM_SPECIES; ++n) {
    const amrex::Real Xface = 0.5*(mole1[n]+mole2[n]);
    const amrex::Real Yface = 0.5*(mass1[n]+mass2[n]);
    const amrex::Real hface = 0.5*(hi1[n]+hi2[n]);
    const amrex::Real dXdy = (mole1[n]-mole2[n])*dxinv[1];
    const amrex::Real Vd = -0.5*(coeffs(iv,CRHOD+n)+coeffs(ivm,CRHOD+n))*(dXdy+(Xface-Yface)*dlnp);
    Vc += Vd;
    fy(iv,UFS+n) += Vd;      
    fy(iv,UEDEN) += Vd*hface;
  }
  // Add correction velocity to fluxes so sum(Vd) = 0
  for (int ns = 0; ns < NUM_SPECIES; ++ns) {
    const amrex::Real Yface = 0.5*(mass1[ns]+mass2[ns]);
    const amrex::Real hface = 0.5*(hi1[ns]+hi2[ns]);
    fy(iv,UFS+ns) -= Yface*Vc;
    fy(iv,UEDEN) -= Yface*hface*Vc;
  }
}
#endif

#if (AMREX_SPACEDIM == 3)
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void
cns_diff_z (int i, int j, int k,
            amrex::Array4<amrex::Real const> const& q,
            amrex::Array4<amrex::Real const> const& coeffs,
            amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dxinv,
            amrex::Array4<amrex::Real> const& fz) noexcept
{
  using amrex::Real;

  Real dTdz = (q(i,j,k,QTEMP)-q(i,j,k-1,QTEMP))*dxinv[2];
  Real dudz = (q(i,j,k,QU)-q(i,j,k-1,QU))*dxinv[2];
  Real dvdz = (q(i,j,k,QV)-q(i,j,k-1,QV))*dxinv[2];
  Real dwdz = (q(i,j,k,QW)-q(i,j,k-1,QW))*dxinv[2];
  Real dudx = (q(i+1,j,k,QU)+q(i+1,j,k-1,QU)-q(i-1,j,k,QU)-q(i-1,j,k-1,QU))*(0.25*dxinv[0]);
  Real dwdx = (q(i+1,j,k,QW)+q(i+1,j,k-1,QW)-q(i-1,j,k,QW)-q(i-1,j,k-1,QW))*(0.25*dxinv[0]);
  Real dvdy = (q(i,j+1,k,QV)+q(i,j+1,k-1,QV)-q(i,j-1,k,QV)-q(i,j-1,k-1,QV))*(0.25*dxinv[1]);
  Real dwdy = (q(i,j+1,k,QW)+q(i,j+1,k-1,QW)-q(i,j-1,k,QW)-q(i,j-1,k-1,QW))*(0.25*dxinv[1]);
  Real divu = dudx + dvdy + dwdz;
  Real mu  = 0.5*(coeffs(i,j,k,CMU)+coeffs(i,j,k-1,CMU));
  Real xi  = 0.5*(coeffs(i,j,k,CXI)+coeffs(i,j,k-1,CXI));
  Real lambda = 0.5*(coeffs(i,j,k,CLAM)+coeffs(i,j,k-1,CLAM));
  Real tauxz = mu*(dudz+dwdx);
  Real tauyz = mu*(dvdz+dwdy);
  Real tauzz = mu*(2.0*dwdz-(2.0/3.0)*divu) + xi*divu;

  fz(i,j,k,UMX)   += -tauxz;
  fz(i,j,k,UMY)   += -tauyz;
  fz(i,j,k,UMZ)   += -tauzz;
  fz(i,j,k,UEDEN) += -0.5*( (q(i,j,k,QU)+q(i,j,k-1,QU))*tauxz
                           +(q(i,j,k,QV)+q(i,j,k-1,QV))*tauyz
                           +(q(i,j,k,QW)+q(i,j,k-1,QW))*tauzz )
                     -lambda*dTdz;

  // Species transport
  const amrex::IntVect iv{AMREX_D_DECL(i, j, k)};
  const amrex::IntVect ivm(iv - amrex::IntVect::TheDimensionVector(2));

  // Get massfrac, molefrac, enthalpy
  Real mass1[NUM_SPECIES], mass2[NUM_SPECIES];
  Real mole1[NUM_SPECIES], mole2[NUM_SPECIES];
  Real hi1[NUM_SPECIES], hi2[NUM_SPECIES];
  for (int n = 0; n < NUM_SPECIES; ++n) {
    mass1[n] = q(iv,QFS+n);
    mass2[n] = q(ivm,QFS+n);
  }
  auto eos = pele::physics::PhysicsType::eos();
  eos.Y2X(mass1, mole1);
  eos.Y2X(mass2, mole2);

  // Compute species and enthalpy fluxes for ideal EOS
  // Get species/enthalpy diffusion, compute correction vel
  Real T = q(iv,QTEMP);
  eos.T2Hi(T, hi1);
  T = q(ivm,QTEMP);
  eos.T2Hi(T, hi2);

  Real Vc = 0.0;
  const amrex::Real dpdz = (q(iv,QPRES)-q(ivm,QPRES))*dxinv[2];
  const amrex::Real dlnp = dpdz / (0.5*(q(iv,QPRES)+q(ivm,QPRES)));
  for (int n = 0; n < NUM_SPECIES; ++n) {
    const amrex::Real Xface = 0.5*(mole1[n]+mole2[n]);
    const amrex::Real Yface = 0.5*(mass1[n]+mass2[n]);
    const amrex::Real hface = 0.5*(hi1[n]+hi2[n]);
    const amrex::Real dXdz = (mole1[n]-mole2[n])*dxinv[2];
    const amrex::Real Vd = -0.5*(coeffs(iv,CRHOD+n)+coeffs(ivm,CRHOD+n))*(dXdz+(Xface-Yface)*dlnp);
    Vc += Vd;
    fz(iv,UFS+n) += Vd;      
    fz(iv,UEDEN) += Vd*hface;
  }
  // Add correction velocity to fluxes so sum(Vd) = 0
  for (int ns = 0; ns < NUM_SPECIES; ++ns) {
    const amrex::Real Yface = 0.5*(mass1[ns]+mass2[ns]);
    const amrex::Real hface = 0.5*(hi1[ns]+hi2[ns]);
    fz(iv,UFS+ns) -= Yface*Vc;
    fz(iv,UEDEN) -= Yface*hface*Vc;
  }
}
#endif

/**
 * @brief Compute diffusion fluxes.
 * 
 * @param i,j,k  x, y, z index. 
 * @param dir    direction, 0:x, 1:y, 2:z.
 * @param q      primitive variables.
 * @param coeffs transport coefficients.
 * @param dxinv  1/dx.
 * @param[out] flx  output diffusion fluxes.
 */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
cns_diff (int i, int j, int k, int dir,
          amrex::Array4<const amrex::Real> const& q,
          amrex::Array4<const amrex::Real> const& coeffs,
          amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dxinv,
          amrex::Array4<amrex::Real> const& flx) noexcept
{
  if (dir == 0) {
    cns_diff_x(i, j, k, q, coeffs, dxinv, flx);
  }
#if (AMREX_SPACEDIM >= 2)
  else if (dir == 1) {
    cns_diff_y(i, j, k, q, coeffs, dxinv, flx);
  }
#endif
#if (AMREX_SPACEDIM == 3)
  else {
    cns_diff_z(i, j, k, q, coeffs, dxinv, flx);
  }
#endif
}

#endif