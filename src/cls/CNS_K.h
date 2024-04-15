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
