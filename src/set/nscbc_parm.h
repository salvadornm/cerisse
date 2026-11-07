#ifndef NSCBC_PARM_H_
#define NSCBC_PARM_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuContainers.H>

namespace nscbc {

enum InflowTargetMode : int {
    TARGET_VELOCITY  = 0,
    TARGET_MASS_FLUX = 1
};

// Default parameter object NSCBC (defined once per boundary).
struct NSCBCParm
{
    amrex::Real Lchar   = 1.0;
    amrex::Real Mmax    = 0.1;
    amrex::Real Ptarget = 101325.0;
    amrex::Real sigma   = 0.28;

    amrex::Real utarget = 0.0;
    amrex::Real vtarget = 0.0;
    amrex::Real wtarget = 0.0;
    amrex::Real Ttarget = 300.0;
    amrex::Real eta     = 1.0;

    // Inflow mass-fraction targets. The default is the first species only.
    amrex::GpuArray<amrex::Real, NUM_SPECIES> Ytarget = {amrex::Real(1.0)};

    bool use_transverse = false;	

    amrex::Real beta_transverse = 1.0;

    // Inflow target selection
    int inflow_target = TARGET_VELOCITY;

    // Positive into the computational domain [kg/(m^2 s)]
    amrex::Real mass_flux_target = 0.0;
};

} // namespace nscbc

#endif
