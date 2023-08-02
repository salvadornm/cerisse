#ifndef CNS_DERIVE_H
#define CNS_DERIVE_H

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>

#include "CNS.H"

void cns_dertemp(const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp,
                 int ncomp, const amrex::FArrayBox& datafab,
                 const amrex::Geometry& geomdata, amrex::Real time, const int* bcrec,
                 int level);

void cns_derpres(const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp,
                 int ncomp, const amrex::FArrayBox& datafab,
                 const amrex::Geometry& geomdata, amrex::Real time, const int* bcrec,
                 int level);

void cns_dereint(const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp,
                 int ncomp, const amrex::FArrayBox& datafab,
                 const amrex::Geometry& geomdata, amrex::Real time, const int* bcrec,
                 int level);

void cns_dervel(const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int ncomp,
                const amrex::FArrayBox& datafab, const amrex::Geometry& geomdata,
                amrex::Real time, const int* bcrec, int level);

void cns_dermagvort(const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp,
                    int ncomp, const amrex::FArrayBox& datafab,
                    const amrex::Geometry& geomdata, amrex::Real time,
                    const int* bcrec, int level);

void cns_dermachnumber(const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp,
                       int ncomp, const amrex::FArrayBox& datafab,
                       const amrex::Geometry& geomdata, amrex::Real time,
                       const int* bcrec, int level);

void cns_derdivu(const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp,
                 int ncomp, const amrex::FArrayBox& datafab,
                 const amrex::Geometry& geomdata, amrex::Real time, const int* bcrec,
                 int level);

void cns_derdivrho(const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp,
                   int ncomp, const amrex::FArrayBox& datafab,
                   const amrex::Geometry& geomdata, amrex::Real time,
                   const int* bcrec, int level);

// void cns_derextsrc (
//   const amrex::Box& bx,
//   amrex::FArrayBox& derfab,
//   int dcomp,
//   int ncomp,
//   const amrex::FArrayBox& datafab,
//   const amrex::Geometry& geomdata,
//   amrex::Real time,
//   const int* bcrec,
//   int level);

void cns_dercp(const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int ncomp,
               const amrex::FArrayBox& datafab, const amrex::Geometry& geomdata,
               amrex::Real time, const int* bcrec, int level);

void cns_dercv(const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int ncomp,
               const amrex::FArrayBox& datafab, const amrex::Geometry& geomdata,
               amrex::Real time, const int* bcrec, int level);

void cns_dermassfrac(const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp,
                     int ncomp, const amrex::FArrayBox& datafab,
                     const amrex::Geometry& geomdata, amrex::Real time,
                     const int* bcrec, int level);

void cns_dermolefrac(const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp,
                     int ncomp, const amrex::FArrayBox& datafab,
                     const amrex::Geometry& geomdata, amrex::Real time,
                     const int* bcrec, int level);

#if NUM_FIELD > 0
void cns_dervaru(const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp,
                 int ncomp, const amrex::FArrayBox& datafab,
                 const amrex::Geometry& geomdata, amrex::Real time, const int* bcrec,
                 int level);

void cns_dertke(const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int ncomp,
                const amrex::FArrayBox& datafab, const amrex::Geometry& geomdata,
                amrex::Real time, const int* bcrec, int level);

void cns_dervary(const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp,
                 int ncomp, const amrex::FArrayBox& datafab,
                 const amrex::Geometry& geomdata, amrex::Real time, const int* bcrec,
                 int level);
#endif

void cns_dervelgrad(const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp,
                    int ncomp, const amrex::FArrayBox& datafab,
                    const amrex::Geometry& geomdata, amrex::Real time,
                    const int* bcrec, int level);

// void cns_dermut(const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int
// ncomp,
//                 const amrex::FArrayBox& datafab, const amrex::Geometry& geomdata,
//                 amrex::Real time, const int* bcrec, int level);

void cns_derturbvisc(const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp,
                     int ncomp, const amrex::FArrayBox& datafab,
                     const amrex::Geometry& geomdata, amrex::Real time,
                     const int* bcrec, int level);

#if CNS_USE_EB
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void get_idx(const int i, const int dir,
                                                 const bool all_regular,
                                                 const amrex::EBCellFlag& flag,
                                                 int& im, int& ip)
{
  if (all_regular) {
    im = i - 1;
    ip = i + 1;
  } else {
    if (flag.isCovered()) {
      im = i;
      ip = i;
    } else {
      const amrex::IntVect ivd = amrex::IntVect::TheDimensionVector(dir);
      im = i - static_cast<int>(flag.isConnected(-ivd));
      ip = i + static_cast<int>(flag.isConnected(ivd));
    }
  }
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real get_weight(const int im,
                                                           const int ip)
{
  const int diff = ip - im;
  return diff == 0 ? 0.0 : (diff == 1 ? 1.0 : 0.5);
}
#endif

#endif