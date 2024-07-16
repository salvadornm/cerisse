#ifndef CNS_PROB_PARM_H
#define CNS_PROB_PARM_H

#include <AMReX_GpuMemory.H>
#include <AMReX_REAL.H>

struct ProbParm
{
  amrex::Real reynolds = 100.0; 
  amrex::Real mach = 0.5;
  amrex::Real prandtl = 0.7;
  amrex::Real p0 = 1.0e5;
  amrex::Real T0 = 3.0e-2;

  amrex::Real L; // length scale (= L_x)  
  amrex::Real rho0;
  amrex::Real u0;

  // MASA parameters
  amrex::Real rho_x_fact = 0.1;
  amrex::Real rho_y_fact = 0.15;
  amrex::Real rho_z_fact = 0.0;
  amrex::Real u_x_fact = 0.06;
  amrex::Real u_y_fact = -0.17;
  amrex::Real u_z_fact = 0.0;
  amrex::Real v_0_fact = 1.3;
  amrex::Real v_x_fact = -0.3;
  amrex::Real v_y_fact = 0.06;
  amrex::Real v_z_fact = 0.0;
  amrex::Real w_0_fact = 0.0;
  amrex::Real w_x_fact = 0.0;
  amrex::Real w_y_fact = 0.0;
  amrex::Real w_z_fact = 0.0;
  amrex::Real p_x_fact = -0.3;
  amrex::Real p_y_fact = 0.2;
  amrex::Real p_z_fact = 0.0;
  amrex::Real a_rhox = 2.0;
  amrex::Real a_rhoy = 4.0;
  amrex::Real a_rhoz = 0.0;
  amrex::Real a_ux = 2.0;
  amrex::Real a_uy = 4.0;
  amrex::Real a_uz = 0.0;
  amrex::Real a_vx = 4.0;
  amrex::Real a_vy = 2.0;
  amrex::Real a_vz = 0.0;
  amrex::Real a_wx = 0.0;
  amrex::Real a_wy = 0.0;
  amrex::Real a_wz = 0.0;
  amrex::Real a_px = 2.0;
  amrex::Real a_py = 4.0;
  amrex::Real a_pz = 0.0;
};

#endif
