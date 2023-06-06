#include <AMReX_Config.H>
#include <AMReX_REAL.H>
#include <parm.H>

#include "diffusion_eb.H"
#include "hyperbolics.H"

// This version works for adiabatic or isothermal wall
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
compute_hyp_wallflux (const amrex::Real rho,
                      AMREX_D_DECL(const amrex::Real u, const amrex::Real v, const amrex::Real w),
                      const amrex::Real p,   const amrex::Real sp[NUM_SPECIES],
                      const amrex::Real axm, const amrex::Real axp,
                      const amrex::Real aym, const amrex::Real ayp,
#if (AMREX_SPACEDIM == 3)
                      const amrex::Real azm, const amrex::Real azp,
#endif
                      amrex::GpuArray<amrex::Real,NVAR>& fluxw, Parm const& parm, 
                      bool eb_isothermal, amrex::Real wall_temp) noexcept
{
  using amrex::Real;
  
  Real apnorm = std::sqrt( AMREX_D_TERM((axm-axp)*(axm-axp), + (aym-ayp)*(aym-ayp), + (azm-azp)*(azm-azp)) );
  if (apnorm <= 0.0) amrex::Abort("compute_hyp_wallflux: we are in trouble.");
  
  Real apnorminv = 1./apnorm;
  Real un = AMREX_D_TERM(u*(axm-axp)*apnorminv, + v*(aym-ayp)*apnorminv, + w*(azm-azp)*apnorminv);
  Real AMREX_D_DECL(momfluxn, momfluxt1, momfluxt2);
  
  if (eb_isothermal) {
    // Fluid states
    Real y[NUM_SPECIES];
    for (int n = 0; n < NUM_SPECIES; ++n) y[n] = sp[n];
  
    auto eos = pele::physics::PhysicsType::eos();
    Real e, cs, T;
    eos.RYP2E(rho, y, p, e);
    e += 0.5*(AMREX_D_TERM(u*u, +v*v, +w*w));
    eos.RPY2Cs(rho, p, y, cs);
    eos.RYP2T(rho, p, y, T);

    // Wall states (extrapolate p and y, use Tw to calculate rhow, ew, csw)
    Real rhow, ew, csw;
    eos.PYT2RE(p, y, 2*wall_temp-T, rhow, ew);
    ew += 0.5*(AMREX_D_TERM(u*u, +v*v, +w*w));
    eos.RPY2Cs(rhow, p, y, csw);

    hllc(rho,  AMREX_D_DECL( un, 0., 0.), p, cs,  e,  y,
         rhow, AMREX_D_DECL(-un, 0., 0.), p, csw, ew, y,
         fluxw[URHO], AMREX_D_DECL(momfluxn, momfluxt1, momfluxt2), 
         fluxw[UEDEN], &fluxw[UFS], parm);

    // amrex::Print() << fluxw[URHO] << " " << momfluxn << "(" << p << ") " << momfluxt1 
    //                << " " << fluxw[UEDEN] << " " << fluxw[UFS] << " " << fluxw[UFS+1] 
    //                << " " << rho << " " << rhow << " " << cs << " " << csw << std::endl;
  } else {
    // Adaibatic wall
    fluxw[URHO] = 0.0;
    momfluxn = p;
    fluxw[UEDEN] = 0.0;
    for (int n = 0; n < NUM_SPECIES; ++n) fluxw[UFS+n] = 0.0;
  }

  AMREX_D_TERM(fluxw[UMX] = (axm-axp) * momfluxn;,
               fluxw[UMY] = (aym-ayp) * momfluxn;,
               fluxw[UMZ] = (azm-azp) * momfluxn;);
}

// // TODO: This version accepts problem specific BC at EB (e.g. jets/isothermal/constant heat flux)
// AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
// compute_hyp_wallflux (int i, int j, int k, amrex::GeometryData const& geomdata, const amrex::Real t,
//                       const amrex::Real rho,
//                       AMREX_D_DECL(const amrex::Real u, const amrex::Real v, const amrex::Real w),
//                       const amrex::Real p,   amrex::Real sp[NUM_SPECIES],
//                       const amrex::Real axm, const amrex::Real axp,
//                       const amrex::Real aym, const amrex::Real ayp,
// #if (AMREX_SPACEDIM == 3)
//                       const amrex::Real azm, const amrex::Real azp,
// #endif
//                       amrex::GpuArray<amrex::Real,NVAR>& fluxw, Parm const& parm,
//                       bool use_prob_eb_state=false,
//                       bool eb_isothermal=false, amrex::Real wall_temp=0) noexcept
// {
//   using amrex::Real;
  
//   Real apnorm = std::sqrt( AMREX_D_TERM((axm-axp)*(axm-axp), + (aym-ayp)*(aym-ayp), + (azm-azp)*(azm-azp)) );
//   if (apnorm <= 0.0) amrex::Abort("compute_hyp_wallflux: we are in trouble.");
  
//   Real apnorminv = 1./apnorm;
//   Real AMREX_D_DECL(momfluxn, momfluxt1, momfluxt2);

//   if (use_prob_eb_state) 
//   {
//     // Fluid states
//     Real un = AMREX_D_TERM(u*(axm-axp)*apnorminv, + v*(aym-ayp)*apnorminv, + w*(azm-azp)*apnorminv);   
   
//     auto eos = pele::physics::PhysicsType::eos();
//     Real e, cs;
//     eos.RYP2E(rho, sp, p, e);
//     e += 0.5*(AMREX_D_TERM(u*u, +v*v, +w*w));
//     eos.RPY2Cs(rho, p, sp, cs);

//     // Wall states
//     Real x[3];
//     const amrex::Real* prob_lo = geomdata.ProbLo();
//     const amrex::Real* dx = geomdata.CellSize();
//     x[0] = prob_lo[0] + (i + 0.5) * dx[0];
//     x[1] = prob_lo[1] + (j + 0.5) * dx[1];
//     x[2] = prob_lo[2] + (k + 0.5) * dx[2];

//     Real rhow, unw, pw, yw[NUM_SPECIES], csw, ew;
//     amrex::Abort("NOT YET IMPLEMENTED");
//     // prob_eb_state(x, t, rho, un, p, sp, e, rhow, unw, pw, yw, ew);

//     hllc(rho,  AMREX_D_DECL(un,  0., 0.), p,  cs,  e,  sp,
//          rhow, AMREX_D_DECL(unw, 0., 0.), pw, csw, ew, yw,
//          fluxw[URHO], AMREX_D_DECL(momfluxn, momfluxt1, momfluxt2), 
//          fluxw[UEDEN], &fluxw[UFS], parm);
//   } 
//   else if (eb_isothermal) 
//   {
//     // Fluid states    
//     Real un = AMREX_D_TERM(u*(axm-axp)*apnorminv, + v*(aym-ayp)*apnorminv, + w*(azm-azp)*apnorminv);
   
//     auto eos = pele::physics::PhysicsType::eos();
//     Real e, cs;
//     eos.RYP2E(rho, sp, p, e);
//     e += 0.5*(AMREX_D_TERM(u*u, +v*v, +w*w));
//     eos.RPY2Cs(rho, p, sp, cs);

//     // Wall states
//     Real rhow, ew, csw;
//     eos.PYT2RE(p, sp, wall_temp, rhow, ew);
//     ew += 0.5*(AMREX_D_TERM(u*u, +v*v, +w*w));
//     eos.RPY2Cs(rhow, p, sp, csw);

//     hllc(rho,  AMREX_D_DECL(un, 0., 0.), p,  cs,  e, sp,
//          rhow, AMREX_D_DECL(un, 0., 0.), p, csw, ew, sp,
//          fluxw[URHO], AMREX_D_DECL(momfluxn, momfluxt1, momfluxt2), 
//          fluxw[UEDEN], &fluxw[UFS], parm);
//   } 
//   else 
//   {
//     fluxw[URHO] = 0.0;
//     momfluxn = p;
//     fluxw[UEDEN] = 0.0;
//     for (int n = 0; n < NUM_SPECIES; ++n) fluxw[UFS+n] = 0.0;
//   }  

//   AMREX_D_TERM(fluxw[UMX] = (axm-axp) * momfluxn;,
//                fluxw[UMY] = (aym-ayp) * momfluxn;,
//                fluxw[UMZ] = (azm-azp) * momfluxn;);
// }

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void 
eb_compute_div (int i, int j, int k, 
                IntVect const& blo, IntVect const& bhi,
                Array4<const Real> const& q, 
                Array4<Real> const& divu,
                AMREX_D_DECL(Array4<Real const> const& u,
                              Array4<Real const> const& v,
                              Array4<Real const> const& w),
                AMREX_D_DECL(Array4<Real> const& fx,
                              Array4<Real> const& fy,
                              Array4<Real> const& fz),
                Array4<EBCellFlag const> const& flag,
                Array4<Real const> const& vfrc,
                Array4<Real const> const& bcent,
                Array4<Real const> const& coefs,
                AMREX_D_DECL(Array4<Real const> const& apx,
                              Array4<Real const> const& apy,
                              Array4<Real const> const& apz),
                AMREX_D_DECL(Array4<Real const> const& fcx,
                              Array4<Real const> const& fcy,
                              Array4<Real const> const& fcz),
                GpuArray<Real,AMREX_SPACEDIM> const& dxinv, Parm const& parm,
                bool do_visc,
                bool eb_no_slip, bool eb_isothermal, const Real eb_wall_temp)
{
  AMREX_D_TERM(bool x_high = (i == bhi[0]);,
               bool y_high = (j == bhi[1]);,
               bool z_high = (k == bhi[2]));
  bool valid_cell = AMREX_D_TERM( (blo[0] <= i) && (i <= bhi[0]),
                               && (blo[1] <= j) && (j <= bhi[1]),
                               && (blo[2] <= k) && (k <= bhi[2]) );

  if (flag(i,j,k).isCovered()) {
    for (int n = 0; n < NVAR; ++n) {
      divu(i,j,k,n) = 0.0;
      if (valid_cell) {
        AMREX_D_TERM(
        fx(i,j,k,n) = 0.; 
        if (x_high) {fx(i+1,j,k,n) = 0.;}
        ,
        fy(i,j,k,n) = 0.;            
        if (y_high) {fy(i,j+1,k,n) = 0.;}
        ,
        fz(i,j,k,n) = 0.;
        if (z_high) {fz(i,j,k+1,n) = 0.;}
        );
      }
    }
  }
  else if (flag(i,j,k).isRegular()) {
    for (int n = 0; n < NVAR; ++n) {
      divu(i,j,k,n) = AMREX_D_TERM(dxinv[0] * (u(i,j,k,n)-u(i+1,j,k,n)) ,
                                 + dxinv[1] * (v(i,j,k,n)-v(i,j+1,k,n)) ,
                                 + dxinv[2] * (w(i,j,k,n)-w(i,j,k+1,n)););
      if (valid_cell) {
        AMREX_D_TERM(
        fx(i,j,k,n) = u(i,j,k,n);
        if (x_high) {fx(i+1,j,k,n) = u(i+1,j,k,n);}
        ,
        fy(i,j,k,n) = v(i,j,k,n);
        if (y_high) {fy(i,j+1,k,n) = v(i,j+1,k,n);}
        ,
        fz(i,j,k,n) = w(i,j,k,n);
        if (z_high) {fz(i,j,k+1,n) = w(i,j,k+1,n);}
        );
      }
    }
  }
  else { // cut-cell
    for (int n = 0; n < NVAR; ++n) {
      // xlo
      Real fxm = u(i,j,k,n);
      if (apx(i,j,k) != 0.0 && apx(i,j,k) != 1.0) {
        AMREX_D_TERM(,
          int jj = j + static_cast<int>(amrex::Math::copysign(1.0_rt, fcx(i,j,k,0)));
          Real fracy = flag(i,j,k).isConnected(0,jj-j,0) ? amrex::Math::abs(fcx(i,j,k,0)) : 0.0_rt;
          ,
          int kk = k + static_cast<int>(amrex::Math::copysign(1.0_rt, fcx(i,j,k,1)));
          Real fracz = flag(i,j,k).isConnected(0,0,kk-k) ? amrex::Math::abs(fcx(i,j,k,1)) : 0.0_rt;
        );
#if AMREX_SPACEDIM == 2
        fxm = (1.0-fracy)*fxm + fracy*u(i,jj,k ,n);
#elif AMREX_SPACEDIM == 3
        fxm = (1.0-fracy)*(1.0-fracz)*fxm 
            + fracy *(1.0-fracz)*u(i,jj,k ,n)
            + fracz *(1.0-fracy)*u(i,j ,kk,n) 
            + fracy * fracz *u(i,jj,kk,n);
#endif
      }
      if (valid_cell) {
        fx(i,j,k,n) = fxm;
      }

      // xhi
      Real fxp = u(i+1,j,k,n);
      if (apx(i+1,j,k) != 0.0 && apx(i+1,j,k) != 1.0) {
        AMREX_D_TERM(,
          int jj = j + static_cast<int>(amrex::Math::copysign(1.0_rt,fcx(i+1,j,k,0)));
          Real fracy = flag(i+1,j,k).isConnected(0,jj-j,0) ? amrex::Math::abs(fcx(i+1,j,k,0)) : 0.0_rt;
          ,
          int kk = k + static_cast<int>(amrex::Math::copysign(1.0_rt,fcx(i+1,j,k,1)));
          Real fracz = flag(i+1,j,k).isConnected(0,0,kk-k) ? amrex::Math::abs(fcx(i+1,j,k,1)) : 0.0_rt;
        );
#if AMREX_SPACEDIM == 2
        fxp = (1.0-fracy)*fxp + fracy *u(i+1,jj,k,n);
#elif AMREX_SPACEDIM == 3
        fxp = (1.0-fracy)*(1.0-fracz)*fxp
                + fracy *(1.0-fracz)*u(i+1,jj,k ,n)
                + fracz *(1.0-fracy)*u(i+1,j ,kk,n)
                + fracy *     fracz *u(i+1,jj,kk,n);
#endif
      }
      if (valid_cell && x_high) {
        fx(i+1,j,k,n) = fxp;
      }

      // ylo
      Real fym = v(i,j,k,n);
      if (apy(i,j,k) != 0.0 && apy(i,j,k) != 1.0) {
        AMREX_D_TERM(,
          int ii = i + static_cast<int>(amrex::Math::copysign(1.0_rt,fcy(i,j,k,0)));
          Real fracx = flag(i,j,k).isConnected(ii-i,0,0) ? amrex::Math::abs(fcy(i,j,k,0)) : 0.0_rt;
          ,
          int kk = k + static_cast<int>(amrex::Math::copysign(1.0_rt,fcy(i,j,k,1)));
          Real fracz = flag(i,j,k).isConnected(0,0,kk-k) ? amrex::Math::abs(fcy(i,j,k,1)) : 0.0_rt;
        );
#if (AMREX_SPACEDIM == 2)
        fym = (1.0-fracx)*fym +  fracx *v(ii,j,k,n);
#else
        fym = (1.0-fracx)*(1.0-fracz)*fym
            +      fracx *(1.0-fracz)*v(ii,j,k ,n)
            +      fracz *(1.0-fracx)*v(i ,j,kk,n)
            +      fracx *     fracz *v(ii,j,kk,n);
#endif
      }
      if (valid_cell) {
        fy(i,j,k,n) = fym;
      }

      // yhi
      Real fyp = v(i,j+1,k,n);
      if (apy(i,j+1,k) != 0.0 && apy(i,j+1,k) != 1.0) {
        AMREX_D_TERM(,
          int ii = i + static_cast<int>(amrex::Math::copysign(1.0_rt,fcy(i,j+1,k,0)));
          Real fracx = flag(i,j+1,k).isConnected(ii-i,0,0) ? amrex::Math::abs(fcy(i,j+1,k,0)) : 0.0_rt;
          ,
          int kk = k + static_cast<int>(amrex::Math::copysign(1.0_rt,fcy(i,j+1,k,1)));
          Real fracz = flag(i,j+1,k).isConnected(0,0,kk-k) ? amrex::Math::abs(fcy(i,j+1,k,1)) : 0.0_rt;
        );
#if (AMREX_SPACEDIM == 2)
        fyp = (1.0-fracx)*fyp + fracx *v(ii,j+1,k,n);
#else
        fyp = (1.0-fracx)*(1.0-fracz)*fyp
            +      fracx *(1.0-fracz)*v(ii,j+1,k ,n)
            +      fracz *(1.0-fracx)*v(i ,j+1,kk,n)
            +      fracx *     fracz *v(ii,j+1,kk,n);
#endif
      }
      if (valid_cell && y_high) {
        fy(i,j+1,k,n) = fyp;
      }

#if (AMREX_SPACEDIM == 3)
      // zlo
      Real fzm = w(i,j,k,n);
      if (apz(i,j,k) != 0.0 && apz(i,j,k) != 1.0) {
        int ii = i + static_cast<int>(amrex::Math::copysign(1.0_rt,fcz(i,j,k,0)));
        int jj = j + static_cast<int>(amrex::Math::copysign(1.0_rt,fcz(i,j,k,1)));
        Real fracx = flag(i,j,k).isConnected(ii-i,0,0) ? amrex::Math::abs(fcz(i,j,k,0)) : 0.0_rt;
        Real fracy = flag(i,j,k).isConnected(0,jj-j,0) ? amrex::Math::abs(fcz(i,j,k,1)) : 0.0_rt;

        fzm = (1.0-fracx)*(1.0-fracy)*fzm
            +      fracx *(1.0-fracy)*w(ii,j ,k,n)
            +      fracy *(1.0-fracx)*w(i ,jj,k,n)
            +      fracx *     fracy *w(ii,jj,k,n);
      }
      if (valid_cell) {
        fz(i,j,k,n) = fzm;
      }

      // zhi
      Real fzp = w(i,j,k+1,n);
      if (apz(i,j,k+1) != 0.0 && apz(i,j,k+1) != 1.0) {
        int ii = i + static_cast<int>(amrex::Math::copysign(1.0_rt,fcz(i,j,k+1,0)));
        int jj = j + static_cast<int>(amrex::Math::copysign(1.0_rt,fcz(i,j,k+1,1)));
        Real fracx = flag(i,j,k+1).isConnected(ii-i,0,0) ? amrex::Math::abs(fcz(i,j,k+1,0)) : 0.0_rt;
        Real fracy = flag(i,j,k+1).isConnected(0,jj-j,0) ? amrex::Math::abs(fcz(i,j,k+1,1)) : 0.0_rt;
        fzp = (1.0-fracx)*(1.0-fracy)*fzp
            +      fracx *(1.0-fracy)*w(ii,j ,k+1,n)
            +      fracy *(1.0-fracx)*w(i ,jj,k+1,n)
            +      fracx *     fracy *w(ii,jj,k+1,n);
      }
      if (valid_cell && z_high) {
        fz(i,j,k+1,n) = fzp;
      }
#endif

      divu(i,j,k,n) = (1.0/vfrc(i,j,k)) *
          ( AMREX_D_TERM(dxinv[0] * (apx(i+1,j,k)*fxp - apx(i,j,k)*fxm) ,
                       + dxinv[1] * (apy(i,j+1,k)*fyp - apy(i,j,k)*fym) ,
                       + dxinv[2] * (apz(i,j,k+1)*fzp - apz(i,j,k)*fzm)) );

      // Monitor problem cell
      if (isnan(divu(i,j,k,n))) {
        std::cout << "div_eb(" << i << "," << j << "," << k << "," << n << ") is nan \n";
        std::cout << "fxp=" << fxp << " fxm=" << fxm << " fyp=" << fyp << " fym=" << fym << "\n";
        std::cout << "up=" << u(i+1,j,k,n) << " um=" << u(i,j,k,n) << " vp=" << v(i,j+1,k,n) << " vm=" << v(i,j,k,n) << "\n";
      }
    } //end loop n

    // Wall fluxes
    GpuArray<Real,NVAR> flux_hyp_wall = {0.0};
    Real sp[NUM_SPECIES];
    for (int n = 0; n < NUM_SPECIES; ++n) sp[n] = q(i,j,k,QFS+n);
    compute_hyp_wallflux(q(i,j,k,QRHO),AMREX_D_DECL(q(i,j,k,QU),q(i,j,k,QV),q(i,j,k,QW)),
                         q(i,j,k,QPRES),sp,
                         apx(i,j,k),apx(i+1,j,k),apy(i,j,k),apy(i,j+1,k),
#if (AMREX_SPACEDIM == 3)
                          apz(i,j,k),apz(i,j,k+1),
#endif
                          flux_hyp_wall,parm, eb_isothermal,eb_wall_temp);

    GpuArray<Real,NVAR> flux_diff_wall = {0.0};
    if (do_visc && (eb_no_slip || eb_isothermal)) {
      compute_diff_wallflux(i,j,k,q,coefs,bcent,
                            apx(i,j,k),apx(i+1,j,k),apy(i,j,k),apy(i,j+1,k),
#if (AMREX_SPACEDIM == 3)
                            apz(i,j,k),apz(i,j,k+1),
#endif
                            flux_diff_wall, eb_no_slip,eb_isothermal,eb_wall_temp);
    }

    for (int n = 0; n < NVAR; ++n) {
      // These assume dx == dy == dz
      divu(i,j,k,n) += flux_hyp_wall[n]*dxinv[0]/vfrc(i,j,k);

#if AMREX_SPACEDIM == 2
      divu(i,j,k,n) += flux_diff_wall[n]*dxinv[0]/vfrc(i,j,k);
#elif AMREX_SPACEDIM == 3
      divu(i,j,k,n) += flux_diff_wall[n]*(dxinv[0]*dxinv[0])/vfrc(i,j,k); //because compute_diff_wallflux is nondimensionalised with dx
#endif

      // The operations following this assume we have returned the negative of the divergence of fluxes.
      divu(i,j,k,n) *= -1.0;
      }
  } //end if cut-cell
}
