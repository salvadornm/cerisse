#ifndef CNS_PROB_H_
#define CNS_PROB_H_

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include "CNS_parm.H"
#include "cns_prob_parm.H"
#include <AMReX_PROB_AMR_F.H>
#include <AMReX_ParmParse.H>

AMREX_GPU_DEVICE inline
void prob_initdata (int i, int j, int k, amrex::Array4<amrex::Real> const& state, amrex::GeometryData const& geomdata, Parm const& parm, ProbParm const& prob_parm) {
  using amrex::Real;
  const Real* prob_lo = geomdata.ProbLo();
  const Real* prob_hi = geomdata.ProbHi();
  const Real* dx      = geomdata.CellSize();

  Real x = prob_lo[0] + (i+Real(0.5))*dx[0];
  Real Pt, rhot, uxt;
  if (x < prob_hi[0]/2) {
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

inline
void user_tagging(amrex::TagBoxArray& tags, amrex::MultiFab &sdata, int level) {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
  Real dengrad_threshold= Real(0.2);
  for (amrex::MFIter mfi(tags,amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
  {
    const amrex::Box& bx = mfi.tilebox();
    auto const& tagfab = tags.array(mfi);
    auto const& sdf = sdata.array(mfi); // state data fab (sdf)
    int idx = 0; // density index
    amrex::ParallelFor(bx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept 
    {
      // amrex::Real ax = amrex::Math::abs(sdf(i+1,j,k,idx) - sdf(i,j,k,idx));
      // amrex::Real ay = amrex::Math::abs(sdf(i,j+1,k,idx) - sdf(i,j,k,idx));
      // amrex::Real az = amrex::Math::abs(sdf(i,j,k+1,idx) - sdf(i,j,k,idx));
      // ax = amrex::max(ax,amrex::Math::abs(sdf(i,j,k,idx) - sdf(i-1,j,k,idx)));
      // ay = amrex::max(ay,amrex::Math::abs(sdf(i,j,k,idx) - sdf(i,j-1,k,idx)));
      // az = amrex::max(az,amrex::Math::abs(sdf(i,j,k,idx) - sdf(i,j,k-1,idx)));
      // if (amrex::max(ax,ay,az) >= dengrad_threshold) {
      //     tagfab(i,j,k) = true;}
      amrex::Real drhox = amrex::Math::abs(sdf(i+1,j,k,idx) - sdf(i-1,j,k,idx))/sdf(i,j,k,idx);
      tagfab(i,j,k) = drhox > 0.5f;
    });
  }
}

  void amrex_probinit (const int* /*init*/, const int* /*name*/, const int* /*namelen*/, const amrex_real* /*problo*/, const amrex_real* /*probhi*/)
  {
    amrex::ParmParse pp("prob");
    pp.query("p_l", CNS::h_prob_parm->p_l);
    pp.query("p_r", CNS::h_prob_parm->p_r);
    pp.query("rho_l", CNS::h_prob_parm->rho_l);
    pp.query("rho_r", CNS::h_prob_parm->rho_r);
    pp.query("u_l", CNS::h_prob_parm->u_l);
    pp.query("u_r", CNS::h_prob_parm->u_r);
#if AMREX_USE_GPU
    amrex::Gpu::htod_memcpy(CNS::d_prob_parm, CNS::h_prob_parm, sizeof(ProbParm));
#endif
  }

AMREX_GPU_DEVICE inline
void user_source(int i, int j, int k, const auto& state, const auto& rhs, const ProbParm& lprobparm) {


}
#endif
