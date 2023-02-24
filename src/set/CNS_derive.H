#ifndef CNS_DERIVE_H_
#define CNS_DERIVE_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>

void cns_derpres (const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int ncomp,
                  const amrex::FArrayBox& datafab, const amrex::Geometry& geomdata,
                  amrex::Real time, const int* bcrec, int level);

// void cns_dervel (const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int ncomp,
//                  const amrex::FArrayBox& datafab, const amrex::Geometry& geomdata,
//                  amrex::Real time, const int* bcrec, int level);
                
void pc_dervelx(
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int dcomp,
  int ncomp,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& geomdata,
  amrex::Real time,
  const int* bcrec,
  const int level);

#endif
