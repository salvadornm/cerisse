#ifndef CNS_PROB_PARM_H
#define CNS_PROB_PARM_H

#include <AMReX_GpuMemory.H>
#include <AMReX_REAL.H>

struct ProbParm
{
  amrex::Real rho;
  amrex::Real u;
  amrex::Real ei;
  amrex::GpuArray<amrex::Real, NUM_SPECIES> Y = {0.0};

  amrex::Real rho_j;
  amrex::Real v_j;
  amrex::Real ei_j;
  amrex::GpuArray<amrex::Real, NUM_SPECIES> Y_jet = {0.0};

  amrex::Real cv_Tinf;

  bool spark = false; // add a hot spot in cavity floor (in bcnormal)

  bool record_statistics = false;    // record time avg (in prob_post_timestep)
  bool clean_aux_on_restart = false; // reset time avg (in prob_post_restart)

  bool do_bl = false;              // add BL at inflow (in bcnormal and prob_init)
  bool make_bl_on_restart = false; // set states before jet to be BL inflow condition (in prob_post_restart)

  bool make_init_on_restart = false; // reset states before x_reset to be inflow condition (in prob_post_restart)
  amrex::Real x_reset; 
};

#endif
