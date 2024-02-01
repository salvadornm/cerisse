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

AMREX_GPU_DEVICE inline
void prob_initdata (int i, int j, int k, amrex::Array4<amrex::Real> const& state, amrex::GeometryData const& geomdata, Parm const& parm, ProbParm const& prob_parm) {
  using amrex::Real;
  const Real* prob_lo = geomdata.ProbLo();
  const Real* prob_hi = geomdata.ProbHi();
  const Real* dx      = geomdata.CellSize();

  Real x = prob_lo[0] + (i+Real(0.5))*dx[0];
  Real Pt, rhot, uxt;
  if (x < prob_hi[0]/10) {
      Pt = prob_parm.p_l;
      rhot = prob_parm.rho_l;
      uxt = prob_parm.u_l;
  } else {
      Pt = prob_parm.p_r;
      rhot = prob_parm.rho_r;
      uxt = prob_parm.u_r;
  }
  state(i,j,k,URHO ) = rhot;
  state(i,j,k,UMX  ) = rhot*uxt;
  state(i,j,k,UMY  ) = Real(0.0);
  state(i,j,k,UMZ  ) = Real(0.0);
  Real et = Pt/(parm.eos_gamma-Real(1.0));
  state(i,j,k,UET) = et + Real(0.5)*rhot*uxt*uxt;
}

void user_tagging(amrex::TagBoxArray& tags, amrex::MultiFab &sdata, int level, IBM::IBMultiFab* ibdata= nullptr);

AMREX_GPU_DEVICE inline
void user_source(int i, int j, int k, const auto& state, const auto& rhs, const ProbParm& lprobparm, Parm const& parms, auto const dx) {
};
#endif