#ifndef CNS_DIFFUSION_K_H_
#define CNS_DIFFUSION_K_H_

#include "index_macros.H"
// #include "parm.H"

#include <AMReX_FArrayBox.H>
#include <AMReX_CONSTANTS.H>
#include <cmath>

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
    Real  dudy = (q(i,j+1,k,QU)+q(i-1,j+1,k,QU)-q(i,j-1,k,QU)-q(i-1,j-1,k,QU))*(0.25*dxinv[1]);
    Real  dvdy = (q(i,j+1,k,QV)+q(i-1,j+1,k,QV)-q(i,j-1,k,QV)-q(i-1,j-1,k,QV))*(0.25*dxinv[1]);
#endif
#if (AMREX_SPACEDIM == 3)
    Real  dudz = (q(i,j,k+1,QU)+q(i-1,j,k+1,QU)-q(i,j,k-1,QU)-q(i-1,j,k-1,QU))*(0.25*dxinv[2]);
    Real  dwdz = (q(i,j,k+1,QW)+q(i-1,j,k+1,QW)-q(i,j,k-1,QW)-q(i-1,j,k-1,QW))*(0.250*dxinv[2]);
#endif
    Real  divu = AMREX_D_TERM(dudx, + dvdy, + dwdz);

    Real  etaf = 0.5*(coeffs(i,j,k,CMU)+coeffs(i-1,j,k,CMU));
    Real  xif  = 0.5*(coeffs(i,j,k,CXI)+coeffs(i-1,j,k,CXI));

    AMREX_D_TERM(Real  tauxx = etaf*(2.0*dudx-(2.0/3.0)*divu) + xif*divu; ,
                 Real  tauxy = etaf*(dudy+dvdx); ,
                 Real  tauxz = etaf*(dudz+dwdx););

    AMREX_D_TERM(fx(i,j,k,UMX)   += -tauxx; ,
                 fx(i,j,k,UMY)   += -tauxy; ,
                 fx(i,j,k,UMZ)   += -tauxz;);

    fx(i,j,k,UEDEN) += -0.5*( AMREX_D_TERM( (q(i,j,k,QU)+q(i-1,j,k,QU))*tauxx,
                                           +(q(i,j,k,QV)+q(i-1,j,k,QV))*tauxy,
                                           +(q(i,j,k,QW)+q(i-1,j,k,QW))*tauxz)
                             + (coeffs(i,j,k,CLAM) +coeffs(i-1,j,k,CLAM))*dTdx);
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

    Real  etaf = 0.5*(coeffs(i,j,k,CMU)+coeffs(i,j-1,k,CMU));
    Real  xif  = 0.5*(coeffs(i,j,k,CXI)+coeffs(i,j-1,k,CXI));
    
    AMREX_D_TERM(Real  tauyy = etaf*(2.0*dvdy-(2.0/3.0)*divu) + xif*divu; ,
                 Real  tauxy = etaf*(dudy+dvdx); ,
                 Real  tauyz = etaf*(dwdy+dvdz););

    AMREX_D_TERM(fy(i,j,k,UMX) += -tauxy; ,
                 fy(i,j,k,UMY) += -tauyy; ,
                 fy(i,j,k,UMZ) += -tauyz;);

    fy(i,j,k,UEDEN) += -0.5*( AMREX_D_TERM( (q(i,j,k,QU)+q(i,j-1,k,QU))*tauxy,
                                           +(q(i,j,k,QV)+q(i,j-1,k,QV))*tauyy,
                                           +(q(i,j,k,QW)+q(i,j-1,k,QW))*tauyz)
                             +(coeffs(i,j,k,CLAM) +coeffs(i,j-1,k,CLAM))*dTdy);
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

    Real  dTdz = (q(i,j,k,QTEMP)-q(i,j,k-1,QTEMP))*dxinv[2];
    Real  dudz = (q(i,j,k,QU)-q(i,j,k-1,QU))*dxinv[2];
    Real  dvdz = (q(i,j,k,QV)-q(i,j,k-1,QV))*dxinv[2];
    Real  dwdz = (q(i,j,k,QW)-q(i,j,k-1,QW))*dxinv[2];
    Real  dudx = (q(i+1,j,k,QU)+q(i+1,j,k-1,QU)-q(i-1,j,k,QU)-q(i-1,j,k-1,QU))*(0.25*dxinv[0]);
    Real  dwdx = (q(i+1,j,k,QW)+q(i+1,j,k-1,QW)-q(i-1,j,k,QW)-q(i-1,j,k-1,QW))*(0.25*dxinv[0]);
    Real  dvdy = (q(i,j+1,k,QV)+q(i,j+1,k-1,QV)-q(i,j-1,k,QV)-q(i,j-1,k-1,QV))*(0.25*dxinv[1]);
    Real  dwdy = (q(i,j+1,k,QW)+q(i,j+1,k-1,QW)-q(i,j-1,k,QW)-q(i,j-1,k-1,QW))*(0.25*dxinv[1]);
    Real  divu = dudx + dvdy + dwdz;
    Real  etaf = 0.5*(coeffs(i,j,k,CMU)+coeffs(i,j,k-1,CMU));
    Real  xif  = 0.5*(coeffs(i,j,k,CXI)+coeffs(i,j,k-1,CXI));
    Real  tauxz = etaf*(dudz+dwdx);
    Real  tauyz = etaf*(dvdz+dwdy);
    Real  tauzz = etaf*(2.0*dwdz-(2.0/3.0)*divu) + xif*divu;

    fz(i,j,k,UMX)   += -tauxz;
    fz(i,j,k,UMY)   += -tauyz;
    fz(i,j,k,UMZ)   += -tauzz;
    fz(i,j,k,UEDEN) += -0.5*( (q(i,j,k,QU)+q(i,j,k-1,QU))*tauxz
                             +(q(i,j,k,QV)+q(i,j,k-1,QV))*tauyz
                             +(q(i,j,k,QW)+q(i,j,k-1,QW))*tauzz
                             +(coeffs(i,j,k,CLAM) +coeffs(i,j,k-1,CLAM))*dTdz );
}
#endif

#endif