#ifndef CNS_BCFILL_H_
#define CNS_BCFILL_H_

void cns_bcfill (
    amrex::Box const& bx, 
    amrex::FArrayBox& data,
    const int dcomp, 
    const int numcomp,
    amrex::Geometry const& geom, 
    const amrex::Real time,
    const amrex::Vector<amrex::BCRec>& bcr, 
    const int bcomp,
    const int scomp);

void cns_react_bcfill (
    amrex::Box const& bx, 
    amrex::FArrayBox& data,
    const int dcomp, 
    const int numcomp,
    amrex::Geometry const& geom, 
    const amrex::Real time,
    const amrex::Vector<amrex::BCRec>& bcr, 
    const int bcomp,
    const int scomp);

#endif