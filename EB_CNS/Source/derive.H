#ifndef CNS_DERIVE_H
#define CNS_DERIVE_H

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>

void cns_derpres(const amrex::Box& bx, 
                 amrex::FArrayBox& derfab, 
                 int dcomp, 
                 int ncomp,
                 const amrex::FArrayBox& datafab, 
                 const amrex::Geometry& geomdata,
                 amrex::Real time, 
                 const int* bcrec, 
                 int level);

void cns_dervel(const amrex::Box& bx, 
                amrex::FArrayBox& derfab, 
                int dcomp, 
                int ncomp,
                const amrex::FArrayBox& datafab, 
                const amrex::Geometry& geomdata,
                amrex::Real time, 
                const int* bcrec, 
                int level);

#endif
