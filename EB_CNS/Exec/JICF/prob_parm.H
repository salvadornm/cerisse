#ifndef CNS_PROB_PARM_H
#define CNS_PROB_PARM_H

#include <AMReX_GpuMemory.H>
#include <AMReX_REAL.H>

struct ProbParm
{
  // Inflow
  amrex::Real rho_inf;
  amrex::Real u_inf;
  amrex::Real ei_inf;
  amrex::Real T_inf;
  amrex::GpuArray<amrex::Real, NUM_SPECIES> Y = {0.0};
  amrex::Real theta0 = 0.5; // initial BL thickness [cm]
  amrex::Real A = 0.05; // perturbation amplitude

  // Jet
  amrex::Real rho_j;
  amrex::Real u_j;
  amrex::Real ei_j;
  amrex::Real T_j;
  amrex::GpuArray<amrex::Real, NUM_SPECIES> Y_j = {0.0};
  amrex::Real r_j = 0.5; // jet radius [cm]

  bool do_spark = false;             // add a hot spot energy deposition (in fill_ext_src)
  bool record_statistics = false;    // record time avg (in prob_post_timestep)
  bool clean_aux_on_restart = false; // reset time avg (in prob_post_restart)
};

#endif
