#ifndef CNS_PROB_PARM_H
#define CNS_PROB_PARM_H

#include <AMReX_GpuMemory.H>
#include <AMReX_REAL.H>

using namespace amrex::literals;

struct ProbParm
{
  // Config 3
  amrex::Real p_1 = 1.5; // top right
  amrex::Real u_1 = 0.0;
  amrex::Real v_1 = 0.0;
  amrex::Real rho_1 = 1.5;
  amrex::Real rhoe_1;
  amrex::Real T_1;

  amrex::Real p_2 = 0.3; // top left
  amrex::Real u_2 = 1.206;
  amrex::Real v_2 = 0.0;
  amrex::Real rho_2 = 0.5323;
  amrex::Real rhoe_2;
  amrex::Real T_2;

  amrex::Real p_3 = 0.029; // bottom left
  amrex::Real u_3 = 1.206;
  amrex::Real v_3 = 1.206;
  amrex::Real rho_3 = 0.138;
  amrex::Real rhoe_3;
  amrex::Real T_3;

  amrex::Real p_4 = 0.3; // bottom right
  amrex::Real u_4 = 0.0;
  amrex::Real v_4 = 1.206;
  amrex::Real rho_4 = 0.5323;
  amrex::Real rhoe_4;
  amrex::Real T_4;
};

#endif
