#ifndef CNS_PROB_H_
#define CNS_PROB_H_

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <CNS_index_macros.H>
#include <CNS_parm.H>
#include <cns_prob_parm.H>
#include <AMReX_TagBox.H>
#include <AMReX_AmrLevel.H>
#include <IBM.H>

AMREX_GPU_DEVICE
void
cns_initdata (int i, int j, int k, amrex::Array4<amrex::Real> const& state,
              amrex::GeometryData const& geomdata, Parm const& parm, ProbParm const& prob_parm);

void tagging(amrex::TagBoxArray& tags, amrex::MultiFab &sdata, int level, IBM::IBMultiFab* ibdata= nullptr);

#endif