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

    AMREX_D_TERM(Real tauxx = mu*(2.0*dudx-(2.0/3.0)*divu) + xi*divu; ,
                 Real tauxy = mu*(dudy+dvdx); ,
                 Real tauxz = mu*(dudz+dwdx););

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

    Real  mu = 0.5*(coeffs(i,j,k,CMU)+coeffs(i,j-1,k,CMU));
    Real  xi  = 0.5*(coeffs(i,j,k,CXI)+coeffs(i,j-1,k,CXI));
    
    AMREX_D_TERM(Real tauyy = mu*(2.0*dvdy-(2.0/3.0)*divu) + xi*divu; ,
                 Real tauxy = mu*(dudy+dvdx); ,
                 Real tauyz = mu*(dwdy+dvdz););

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
    Real  mu = 0.5*(coeffs(i,j,k,CMU)+coeffs(i,j,k-1,CMU));
    Real  xi  = 0.5*(coeffs(i,j,k,CXI)+coeffs(i,j,k-1,CXI));
    Real  tauxz = mu*(dudz+dwdx);
    Real  tauyz = mu*(dvdz+dwdy);
    Real  tauzz = mu*(2.0*dwdz-(2.0/3.0)*divu) + xi*divu;

    fz(i,j,k,UMX)   += -tauxz;
    fz(i,j,k,UMY)   += -tauyz;
    fz(i,j,k,UMZ)   += -tauzz;
    fz(i,j,k,UEDEN) += -0.5*( (q(i,j,k,QU)+q(i,j,k-1,QU))*tauxz
                             +(q(i,j,k,QV)+q(i,j,k-1,QV))*tauyz
                             +(q(i,j,k,QW)+q(i,j,k-1,QW))*tauzz
                             +(coeffs(i,j,k,CLAM) +coeffs(i,j,k-1,CLAM))*dTdz );
}
#endif

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void
cns_diff (int i, int j, int k, int dir,
          amrex::Array4<amrex::Real const> const& q,
          amrex::Array4<amrex::Real const> const& coeffs,
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

//     using amrex::Real;

//     const amrex::IntVect iv(AMREX_D_DECL(i, j, k));
//     amrex::IntVect iv_dir[AMREX_SPACEDIM]; 
//     AMREX_D_TERM(iv_dir[0] = amrex::IntVect::TheDimensionVector(0); ,
//                  iv_dir[1] = amrex::IntVect::TheDimensionVector(1); ,
//                  iv_dir[2] = amrex::IntVect::TheDimensionVector(2););
//     int dir1, dir2, QU1, QU2, QU3;
//     if      (dir == 0) { dir1 = 1; dir2 = 2; AMREX_D_TERM(QU1 = QU; QU2 = QV; QU3 = QW;); }
//     else if (dir == 1) { dir1 = 0; dir2 = 2; AMREX_D_TERM(QU1 = QV; QU2 = QU; QU3 = QW;); }
//     else               { dir1 = 0; dir2 = 1; AMREX_D_TERM(QU1 = QW; QU2 = QU; QU3 = QV;); }

//     Real dTdx = (q(iv,QTEMP)-q(iv-iv_dir[dir],QTEMP))*dxinv[dir];
//     AMREX_D_TERM(Real dudx = (q(iv,QU1)-q(iv-iv_dir[dir],QU1))*dxinv[dir];,
//                  Real dvdx = (q(iv,QU2)-q(iv-iv_dir[dir],QU2))*dxinv[dir];,
//                  Real dwdx = (q(iv,QU3)-q(iv-iv_dir[dir],QU3))*dxinv[dir];);
// #if (AMREX_SPACEDIM >= 2)
//     Real dudy = (q(iv+iv_dir[dir1],QU1)+q(iv-iv_dir[dir]+iv_dir[dir1],QU1)
//                 -q(iv-iv_dir[dir1],QU1)-q(iv-iv_dir[dir]-iv_dir[dir1],QU1))*(0.25*dxinv[dir1]);
//     Real dvdy = (q(iv+iv_dir[dir1],QU2)+q(iv-iv_dir[dir]+iv_dir[dir1],QU2)
//                 -q(iv-iv_dir[dir1],QU2)-q(iv-iv_dir[dir]-iv_dir[dir1],QU2))*(0.25*dxinv[dir1]);
// #endif
// #if (AMREX_SPACEDIM == 3)
//     Real dudz = (q(iv+iv_dir[dir2],QU1)+q(iv-iv_dir[dir]+iv_dir[dir2],QU1)
//                 -q(iv-iv_dir[dir2],QU1)-q(iv-iv_dir[dir]-iv_dir[dir2],QU1))*(0.25*dxinv[dir2]);
//     Real dwdz = (q(iv+iv_dir[dir2],QU3)+q(iv-iv_dir[dir]+iv_dir[dir2],QU3)
//                 -q(iv-iv_dir[dir2],QU3)-q(iv-iv_dir[dir]-iv_dir[dir2],QU3))*(0.25*dxinv[dir2]);
// #endif
//     Real divu = AMREX_D_TERM(dudx, +dvdy, +dwdz);

//     Real mu = 0.5*(coeffs(iv,CMU)+coeffs(iv-iv_dir[dir],CMU));
//     Real xi = 0.5*(coeffs(iv,CXI)+coeffs(iv-iv_dir[dir],CXI));

//     AMREX_D_TERM(Real tauxx = mu*(2.0*dudx-(2.0/3.0)*divu) + xi*divu; ,
//                  Real tauxy = mu*(dudy+dvdx); ,
//                  Real tauxz = mu*(dudz+dwdx););

//     AMREX_D_TERM(flx(iv,UMX) += -tauxx; ,
//                  flx(iv,UMY) += -tauxy; ,
//                  flx(iv,UMZ) += -tauxz;);

//     flx(iv,UEDEN) -= 0.5*( AMREX_D_TERM( (q(iv,QU)+q(iv-iv_dir[dir],QU))*tauxx,
//                                         +(q(iv,QV)+q(iv-iv_dir[dir],QV))*tauxy,
//                                         +(q(iv,QW)+q(iv-iv_dir[dir],QW))*tauxz)
//                           + (coeffs(iv,CLAM) +coeffs(iv-iv_dir[dir],CLAM))*dTdx);
}

#endif