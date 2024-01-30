#ifndef CNS_PROB_PARM_H
#define CNS_PROB_PARM_H

#include <AMReX_GpuMemory.H>
#include <AMReX_REAL.H>

struct ProbParm
{
  amrex::Real reynolds = 1600.0; 
  amrex::Real mach = 0.3;
  amrex::Real prandtl = 0.71;
  amrex::Real p0 = 1.013e6;
  amrex::Real T0 = 300.0;
  amrex::Real L = 0.0; // length scale (= L_x)
  
  amrex::Real rho0 = 0.0;
  amrex::Real u0 = 0.0;

  // MASA parameters
  amrex::Real rho_x_fact = 0.1;
  amrex::Real rho_y_fact = 0.15;
  amrex::Real rho_z_fact = -0.05;
  amrex::Real u_x_fact = 0.06;
  amrex::Real u_y_fact = -0.17;
  amrex::Real u_z_fact = 0.0;
  amrex::Real v_0_fact = 1.3;
  amrex::Real v_x_fact = -0.3;
  amrex::Real v_y_fact = 0.05;
  amrex::Real v_z_fact = 0.0;
  amrex::Real w_0_fact = 0.0;
  amrex::Real w_x_fact = 0.0;
  amrex::Real w_y_fact = 0.0;
  amrex::Real w_z_fact = -0.1;
  amrex::Real p_x_fact = -0.3;
  amrex::Real p_y_fact = 0.2;
  amrex::Real p_z_fact = 0.0;
  amrex::Real a_rhox = 2.0;
  amrex::Real a_rhoy = 4.0;
  amrex::Real a_rhoz = 6.0;
  amrex::Real a_ux = 2.0;
  amrex::Real a_uy = 4.0;
  amrex::Real a_uz = 0.0;
  amrex::Real a_vx = 4.0;
  amrex::Real a_vy = 2.0;
  amrex::Real a_vz = 0.0;
  amrex::Real a_wx = 0.0;
  amrex::Real a_wy = 0.0;
  amrex::Real a_wz = 2.0;
  amrex::Real a_px = 6.0;
  amrex::Real a_py = 2.0;
  amrex::Real a_pz = 4.0;
};

#endif
