#ifndef CNS_K_H_
#define CNS_K_H_

#include <AMReX_FArrayBox.H>
#include <prob.h>

#include <cmath>
#include <limits>

using namespace amrex;

inline void derpres(const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                    const FArrayBox& datafab, const Geometry& /*geomdata*/,
                    Real /*time*/, const int* /*bcrec*/, int /*level*/) {
  auto const dat = datafab.const_array();
  auto pfab = derfab.array(dcomp);
  PROB::ProbClosures const* cls = CNS::d_prob_closures;

  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {

    Real Q[cls->NPRIM],U[cls->NCONS];
    for (int n = 0; n < cls->NCONS; ++n) {
      U[n] = dat(i,j,k,n);
    }  
    cls->cons2prims_point(U,Q);
    pfab(i, j, k) = Q[cls->QPRES];    
  });
}

inline void dertemp(const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                    const FArrayBox& datafab, const Geometry& /*geomdata*/,
                    Real /*time*/, const int* /*bcrec*/, int /*level*/) {
  auto const dat = datafab.const_array();
  auto Tfab = derfab.array(dcomp);
  PROB::ProbClosures const* cls = CNS::d_prob_closures;

  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {

    Real Q[cls->NPRIM],U[cls->NCONS];
    for (int n = 0; n < cls->NCONS; ++n) {
      U[n] = dat(i,j,k,n);
    }  
    cls->cons2prims_point(U,Q);
    Tfab(i, j, k) = Q[cls->QT];
  });
}

inline void dervel(const Box& bx, FArrayBox& derfab, int dcomp, int ncomp,
                   const FArrayBox& datfab, const Geometry& /*geomdata*/,
                   Real /*time*/, const int* /*bcrec*/, const int /*level*/) {
  auto const dat = datfab.const_array();
  auto vel = derfab.array();
  PROB::ProbClosures const* cls = CNS::d_prob_closures;

  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    Real rho = 0.0;
    for (int n = 0; n < NUM_SPECIES; ++n) {
      rho += dat(i, j, k, cls->UFS + n);
    }
    Real rhoinv = 1.0_rt / rho;

    AMREX_D_TERM(vel(i, j, k, dcomp) = dat(i, j, k, cls->UMX) * rhoinv;
                 , vel(i, j, k, dcomp + 1) = dat(i, j, k, cls->UMY) * rhoinv;
                 , vel(i, j, k, dcomp + 2) = dat(i, j, k, cls->UMZ) * rhoinv;)
  });
}

inline void derdensity(const Box& bx, FArrayBox& derfab, int dcomp, int ncomp,
                       const FArrayBox& datfab, const Geometry& /*geomdata*/,
                       Real /*time*/, const int* /*bcrec*/,
                       const int /*level*/) {
  auto const dat = datfab.const_array();
  auto den = derfab.array();
  PROB::ProbClosures const* cls = CNS::d_prob_closures;

  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    Real rho = 0.0;
    for (int n = 0; n < NUM_SPECIES; ++n) {
      rho += dat(i, j, k, cls->UFS + n);
    }
    den(i, j, k, dcomp) = rho;
  });
}

inline void derkineticenergy(const Box& bx, FArrayBox& derfab, int dcomp,
                             int ncomp, const FArrayBox& datfab,
                             const Geometry& /*geomdata*/, Real /*time*/,
                             const int* /*bcrec*/, const int /*level*/) {
  auto const dat = datfab.const_array();
  auto ke = derfab.array(dcomp);
  PROB::ProbClosures const* cls = CNS::d_prob_closures;

  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    Real rho = 0.0;
    for (int n = 0; n < NUM_SPECIES; ++n) {
      rho += dat(i, j, k, cls->UFS + n);
    }
    Real mx = dat(i, j, k, cls->UMX);
    Real my = dat(i, j, k, cls->UMY);
    Real mz = dat(i, j, k, cls->UMZ);
    ke(i, j, k) = Real(0.5) * (mx * mx + my * my + mz * mz) / rho;
  });
}

inline void dermagvort(const Box& bx, FArrayBox& derfab, int dcomp, int ncomp,
                       const FArrayBox& datfab, const Geometry& geomdata,
                       Real /*time*/ = -1, const int* /*bcrec*/ = nullptr,
                       const int /*level*/ = -1) {
  auto const dat = datfab.const_array();
  auto vort = derfab.array(dcomp);

  const amrex::Box& gbx = amrex::grow(bx, 1);

  amrex::FArrayBox local(gbx, 3, The_Async_Arena());
  auto larr = local.array();

  PROB::ProbClosures const* cls = CNS::d_prob_closures;

  // Convert momentum to velocity
  amrex::ParallelFor(gbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    Real rho = 0.0;
    for (int n = 0; n < NUM_SPECIES; ++n) {
      rho += dat(i, j, k, cls->UFS + n);
    }
    AMREX_D_TERM(larr(i, j, k, 0) = dat(i, j, k, cls->UMX) / rho;
                 , larr(i, j, k, 1) = dat(i, j, k, cls->UMY) / rho;
                 , larr(i, j, k, 2) = dat(i, j, k, cls->UMZ) / rho;)
  });

  AMREX_D_TERM(const amrex::Real dx = geomdata.CellSize(0);
               , const amrex::Real dy = geomdata.CellSize(1);
               , const amrex::Real dz = geomdata.CellSize(2);)

  // Calculate vorticity
  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    AMREX_D_TERM(vort(i, j, k) = 0.0;
                 , 
                 int im = i - 1; int ip = i + 1; amrex::Real wi = 0.5;
                 int jm = j - 1; int jp = j + 1; amrex::Real wj = 0.5;
                 const amrex::Real vx = wi * (larr(ip, j, k, 1) - larr(im, j, k, 1)) / dx;
                 const amrex::Real uy = wj * (larr(i, jp, k, 0) - larr(i, jm, k, 0)) / dy;
                 const amrex::Real v3 = vx - uy;
                 , 
                 int km = k - 1; int kp = k + 1; amrex::Real wk = 0.5;
                 const amrex::Real wx =  wi * (larr(ip, j, k, 2) - larr(im, j, k, 2)) / dx;
                 const amrex::Real wy = wj * (larr(i, jp, k, 2) - larr(i, jm, k, 2)) / dy;
                 const amrex::Real uz = wk * (larr(i, j, kp, 0) - larr(i, j, km, 0)) / dz;
                 const amrex::Real vz = wk * (larr(i, j, kp, 1) - larr(i, j, km, 1)) / dz;
                 const amrex::Real v1 = wy - vz;
                 const amrex::Real v2 = uz - wx;)
    vort(i, j, k) = std::sqrt(AMREX_D_TERM(0., +v3 * v3, +v1 * v1 + v2 * v2));
  });
}

inline void derenstrophy(const Box& bx, FArrayBox& derfab, int dcomp, int ncomp,
                         const FArrayBox& datfab, const Geometry& geomdata,
                         Real /*time*/, const int* /*bcrec*/,
                         const int /*level*/) {
  auto const dat = datfab.const_array();
  auto ens = derfab.array(dcomp);

  amrex::FArrayBox vortfab(bx, 1, The_Async_Arena());
  dermagvort(bx, vortfab, 0, 1, datfab, geomdata);
  
  PROB::ProbClosures const* cls = CNS::d_prob_closures;

  auto vort = vortfab.const_array();
  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    Real rho = 0.0;
    for (int n = 0; n < NUM_SPECIES; ++n) {
      rho += dat(i, j, k, cls->UFS + n);
    }
    ens(i, j, k) = 0.5 * rho * vort(i, j, k) * vort(i, j, k);
  });
}
#endif
