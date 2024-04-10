#ifndef CNS_K_H_
#define CNS_K_H_

#include <AMReX_FArrayBox.H>
#include <prob.h>

#include <cmath>
#include <limits>

using namespace amrex;

// TODO: Pass primsmf here to avoid computing primitive variables again.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real
cns_estdt(Box const& bx, Array4<Real const> const& state,
          GpuArray<Real, AMREX_SPACEDIM> const& dx,
          const PROB::ProbClosures& cls) noexcept {
  const auto lo = lbound(bx);
  const auto hi = ubound(bx);
#if !defined(__CUDACC__) || (__CUDACC_VER_MAJOR__ != 9) || \
    (__CUDACC_VER_MINOR__ != 2)
  Real dt = std::numeric_limits<Real>::max();
#else
  Real dt = Real(1.e37);
#endif

  for (int k = lo.z; k <= hi.z; ++k) {
    for (int j = lo.y; j <= hi.y; ++j) {
      for (int i = lo.x; i <= hi.x; ++i) {
        Real rho = 0.0;
        for (int n = 0; n < NUM_SPECIES; ++n) {
          rho += state(i, j, k, cls.UFS + n);
        }
        Real mx = state(i, j, k, cls.UMX);
        Real my = state(i, j, k, cls.UMY);
        Real mz = state(i, j, k, cls.UMZ);
        Real rhoinv = Real(1.0) / max(rho, 1.0e-40);
        Real vx = mx * rhoinv;
        Real vy = my * rhoinv;
        Real vz = mz * rhoinv;
        Real ke = Real(0.5) * rho * (vx * vx + vy * vy + vz * vz);
        Real ei = (state(i, j, k, cls.UET) - ke) * rhoinv;
        Real Y[NUM_SPECIES];
        for (int n = 0; n < NUM_SPECIES; ++n) {
          Y[n] = state(i, j, k, cls.UFS + n) * rhoinv;
          // if (Y[n] < 0.0) {
          //   std::cout << "Y[" << n << "](" << i << ", " << j << ", " << k
          //             << ") = " << Y[n] << std::endl;
          // }
        }
        Real cs;
        cls.RYE2Cs(rho, Y, ei, cs);
        Real dtx = dx[0] / (Math::abs(vx) + cs);
        Real dty = dx[1] / (Math::abs(vy) + cs);
        Real dtz = dx[2] / (Math::abs(vz) + cs);
        dt = amrex::min(dt, AMREX_D_DECL(dtx, dty, dtz));
      }
    }
  }

  return dt;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void pointCFL(
    int i, int j, int k, Array2D<Real, 0, 2, 0, 2>& array,
    const Array4<const Real>& prims, const PROB::ProbClosures& cls,
    const auto dx, Real dt) noexcept {
  Real sos = prims(i, j, k, cls.QC);
  //
  for (int idir = 0; idir < AMREX_SPACEDIM; idir++) {
    array(idir, 0) =
        max(array(idir, 0), std::abs(prims(i, j, k, cls.QU + idir) +
                                     sos));  // max(|ui + c|) in +ith direction
    array(idir, 1) =
        max(array(idir, 1), std::abs(prims(i, j, k, cls.QU + idir) -
                                     sos));  // max(|ui - c|) in -ith direction
    array(idir, 2) = max(array(idir, 0), array(idir, 1)) * dt /
                     dx[idir];  // CFL max in ith direction
    // Print() << "inside " << array(idir,0) << " " << array(idir,1) << " " <<
    // array(idir,2) << std::endl;
  }
}

// AMREX_GPU_DEVICE AMREX_FORCE_INLINE void cns_compute_temperature(
//     int i, int j, int k, Array4<Real> const& u,
//     const PROB::ProbClosures& cls) noexcept {
//   Real rhoinv = Real(1.0) / u(i, j, k, cls.URHO);
//   Real mx = u(i, j, k, cls.UMX);
//   Real my = u(i, j, k, cls.UMY);
//   Real mz = u(i, j, k, cls.UMZ);
//   Real rhoeint =
//       u(i, j, k, cls.UET) - Real(0.5) * rhoinv * (mx * mx + my * my + mz * mz);
//   // u(i,j,k,UTEMP) = rhoinv * rhoeint * (Real(1.0)/parm.cv);
// }

inline void derpres(const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                    const FArrayBox& datafab, const Geometry& /*geomdata*/,
                    Real /*time*/, const int* /*bcrec*/, int /*level*/) {
  auto const dat = datafab.const_array();
  auto pfab = derfab.array(dcomp);
  PROB::ProbClosures const* cls = CNS::d_prob_closures;

  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    // cons2P
    Real rho = 0.0;
    for (int n = 0; n < NUM_SPECIES; ++n) {
      rho += dat(i, j, k, cls->UFS + n);
    }
    Real rhoinv = 1.0_rt / rho;
    Real mx = dat(i, j, k, cls->UMX);
    Real my = dat(i, j, k, cls->UMY);
    Real mz = dat(i, j, k, cls->UMZ);
    Real ei = rhoinv * (dat(i, j, k, cls->UET) -
                        Real(0.5) * rhoinv * (mx * mx + my * my + mz * mz));
    Real massfrac[NUM_SPECIES];
    for (int n = 0; n < NUM_SPECIES; ++n) {
      massfrac[n] = dat(i, j, k, cls->UFS + n) * rhoinv;
    }

    Real dummy, pres;
    cls->RYE2TP(rho, massfrac, ei, dummy, pres);
    pfab(i, j, k) = pres;
  });
}

inline void dertemp(const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                    const FArrayBox& datafab, const Geometry& /*geomdata*/,
                    Real /*time*/, const int* /*bcrec*/, int /*level*/) {
  auto const dat = datafab.const_array();
  auto Tfab = derfab.array(dcomp);
  PROB::ProbClosures const* cls = CNS::d_prob_closures;

  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    // cons2T
    Real rho = 0.0;
    for (int n = 0; n < NUM_SPECIES; ++n) {
      rho += dat(i, j, k, cls->UFS + n);
    }    
    Real rhoinv = 1.0_rt / rho;
    Real mx = dat(i, j, k, cls->UMX);
    Real my = dat(i, j, k, cls->UMY);
    Real mz = dat(i, j, k, cls->UMZ);
    Real ei = rhoinv * (dat(i, j, k, cls->UET) -
                        Real(0.5) * rhoinv * (mx * mx + my * my + mz * mz));
    Real massfrac[NUM_SPECIES];
    for (int n = 0; n < NUM_SPECIES; ++n) {
      massfrac[n] = dat(i, j, k, cls->UFS + n) * rhoinv;
    }

    Real temp, dummy;
    cls->RYE2TP(rho, massfrac, ei, temp, dummy);
    Tfab(i, j, k) = temp;
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
#endif
