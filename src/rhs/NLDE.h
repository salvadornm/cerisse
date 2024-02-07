
#include <AMReX_Vector.H>

namespace NLDE {
// vector of baseflow multifabs
inline Vector<MultiFab> Vbaseflow;

// allocate baseflow multifabs
inline void allocVMF(int& nlevs) { Vbaseflow.resize(nlevs); }

// baseflow interpolation
// inline void interpolate_baseflow() {
// }

inline void post_regrid(const int& lev, const BoxArray& grids,
                        const DistributionMapping& dmap, const MFInfo& info,
                        const FabFactory<FArrayBox>& factory) {
  // reallocate MF
  Vbaseflow[lev].clear();
  Vbaseflow[lev].define(grids, dmap, NCONS, 0, info, factory);
  Vbaseflow[lev].setVal(0.0);

  // interpolate_baseflow();
}

// convert primitive variables to conserved variables
inline void cons2prim(const Box& bxg, const Array4<const Real>& u,
                      const Array4<Real>& q, const PROB::ProbClosures& cls) {
  amrex::ParallelFor(bxg, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    Abort("NLDE cons2prim to implement");
  });
}

// computes Euler fluxes of linear disturbances from conserved variable vector
// and maxeigen value
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void eflux_linear(
    int i, int j, int k, auto const& base, auto const& cons, auto const& fx,
    auto const& fy, auto const& fz, auto const& lambda,
    const PROB::ProbClosures& cls) {
  Abort("NLDE Euler fluxes not implemented yet");

  // disturbance fluxes from eqn 15
  // also need baseflow quantities and fluxes
  //   Real rhoinv = Real(1.0)/cons(i,j,k,URHO);
  //   Real momx   = cons(i,j,k,UMX);
  //   Real momy   = cons(i,j,k,UMY);
  //   Real momz   = cons(i,j,k,UMZ);
  //   Real rhoet  = cons(i,j,k,UET);
  //   Real ux     = momx*rhoinv;
  //   Real uy     = momy*rhoinv;
  //   Real uz     = momz*rhoinv;

  //   Real rhoekin = Real(0.5)*rhoinv*(momx*momx + momy*momy + momz*momz);
  //   Real rhoeint = rhoet - rhoekin;
  //   Real P       = (closures.gamma - Real(1.0))*rhoeint;

  //   fx(i,j,k,URHO)  = momx;
  //   fx(i,j,k,UMX)   = momx*ux + P;
  //   fx(i,j,k,UMY)   = momy*ux;
  //   fx(i,j,k,UMZ)   = momz*ux;
  //   fx(i,j,k,UET)   = (rhoet + P)*ux;

  //   fy(i,j,k,URHO)  = momy;
  //   fy(i,j,k,UMX)   = momx*uy;
  //   fy(i,j,k,UMY)   = momy*uy + P;
  //   fy(i,j,k,UMZ)   = momz*uy;
  //   fy(i,j,k,UET)   = (rhoet + P)*uy;

  //   fz(i,j,k,URHO)  = momz;
  //   fz(i,j,k,UMX)   = momx*uz;
  //   fz(i,j,k,UMY)   = momy*uz;
  //   fz(i,j,k,UMZ)   = momz*uz + P;
  //   fz(i,j,k,UET)   = (rhoet + P)*uz;

  //   Real cs=sqrt(closures.gamma*P*rhoinv);
  //   lambda(i,j,k,0) = std::abs(max(ux+cs,ux-cs,ux));
  //   lambda(i,j,k,1) = std::abs(max(uy+cs,uy-cs,uy));
  //   lambda(i,j,k,2) = std::abs(max(uz+cs,uz-cs,uz));
}

// euler flux computation function
AMREX_FORCE_INLINE void eflux(const Box& bxg, const Array4<const Real>& consfab,
                              const Array4<Real>& primfab,
                              const Array4<Real>& pflx,
                              const Array4<Real>& nflx,
                              const PROB::ProbClosures& cls) {
  // loop over all fabs
  // if statement for linear disturance flux calculation (dist_linear)
  // ParallelFor(bxg,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
  //     eflux_linear(i, j, k, basefab, statefab, pflx, lambda, lcls);
  // });

  // ParallelFor(bxg,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
  //     cons2char(i, j, k, statefab, pfabfx, pfabfy, pfabfz, lclosures);
  // });

  // ParallelFor(bxnodal, int(NCONS) , [=] AMREX_GPU_DEVICE (int i, int j, int
  // k, int n) noexcept {
  //       numericalflux_globallaxsplit(i, j, k, n, statefab ,pfabfx, pfabfy,
  //       pfabfz, lambda ,nfabfx, nfabfy, nfabfz);
  // });
}

// similarly, viscous fluxes
AMREX_FORCE_INLINE void vflux(MultiFab& statemf, MultiFab& primsmf,
                              Array<MultiFab, AMREX_SPACEDIM>& numflxmf) {}

}  // namespace NLDE