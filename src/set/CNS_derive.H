#ifndef CNS_DERIVE_H_
#define CNS_DERIVE_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>

void derpres (const amrex::Box& bx, amrex::FArrayBox& derfab, 
                  int dcomp, int ncomp, const amrex::FArrayBox& datafab, const amrex::Geometry& geomdata, amrex::Real time, const int* bcrec, int level);

void dervelx (const amrex::Box& bx, amrex::FArrayBox& derfab, 
              int dcomp, int ncomp, const amrex::FArrayBox& datafab, const amrex::Geometry& geomdata, amrex::Real time, const int* bcrec, int level);

void dervely (const amrex::Box& bx, amrex::FArrayBox& derfab, 
              int dcomp, int ncomp, const amrex::FArrayBox& datafab, const amrex::Geometry& geomdata, amrex::Real time, const int* bcrec, int level);

void dervelz (const amrex::Box& bx, amrex::FArrayBox& derfab, 
              int dcomp, int ncomp, const amrex::FArrayBox& datafab, const amrex::Geometry& geomdata, amrex::Real time, const int* bcrec, int level);
#endif
