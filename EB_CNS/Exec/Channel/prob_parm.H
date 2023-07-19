#ifndef CNS_PROB_PARM_H
#define CNS_PROB_PARM_H

#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>

using namespace amrex;

struct ProbParm
{
  Real H = 0.6845;     // channel half height [cm]
  Real Re_tau = 220.0; // Re_tau = rho_w*u_tau*H/mu_w
  Real Re_b = 3000.0;  // bulk Re = rho_b*u_b*H/mu_w
  Real M_b = 1.5;      // bulk M = u_b/c_w
  Real u_tau = 3500.0; // friction velocity = sqrt(tau_w/rho_w) [cm/s]
  Real T_w = 500.0;    // wall temperature [K]
  Real mu_w = 2.67e-4; // dynamic viscosity at Tw [g/cm-s]
  
  Real rho_w;
  Real u_b;
  Real rho_b;
  Real f_x;

  GpuArray<Real, NUM_SPECIES> massfrac = {1.0};
};

#endif
