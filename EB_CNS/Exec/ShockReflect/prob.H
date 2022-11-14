#ifndef CNS_PROB_H
#define CNS_PROB_H

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_REAL.H>

#include "CNS.H"
#include "prob_parm.H"

/**
 * \brief Initialise state data.
 * 
 * @param i         x position.
 * @param j         y position.
 * @param k         z position.
 * @param state     state data.
 * @param geomdata  domain geometry data.
 * @param parm      Parm data defined in parm.H.
 * @param prob_parm ProbParm data as defined in prob_parm.H and initialised in amrex_probinit.
 */
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
cns_initdata (int i, int j, int k, 
              amrex::Array4<amrex::Real> const& state,
              amrex::GeometryData const& geomdata,
              Parm const& /*parm*/,
              ProbParm const& prob_parm)
{
    // Geometry
    const amrex::Real* prob_lo = geomdata.ProbLo();
    const amrex::Real* dx = geomdata.CellSize();
    const amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];
    const amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1];

    // Set the states
    if (x < 1./6. + 1./sqrt(3.)*y) {
        state(i, j, k, URHO) = prob_parm.rho_l;
        state(i, j, k, UMX) = prob_parm.rho_l * prob_parm.u_l;
        state(i, j, k, UMY) = prob_parm.rho_l * prob_parm.v_l;
        state(i, j, k, UMZ) = 0.0;
        state(i, j, k, UEDEN) = prob_parm.rhoe_l + 0.5 * prob_parm.rho_l * (prob_parm.u_l * prob_parm.u_l + prob_parm.v_l * prob_parm.v_l);
        state(i, j, k, UEINT) = prob_parm.rhoe_l;
        state(i, j, k, UTEMP) = prob_parm.T_l;
        state(i, j, k, UFS) = prob_parm.rho_l;
    } else {
        state(i, j, k, URHO) = prob_parm.rho_r;
        state(i, j, k, UMX) = 0.0;
        state(i, j, k, UMY) = 0.0;
        state(i, j, k, UMZ) = 0.0;
        state(i, j, k, UEDEN) = prob_parm.rhoe_r;
        state(i, j, k, UEINT) = prob_parm.rhoe_r;
        state(i, j, k, UTEMP) = prob_parm.T_r;
        state(i, j, k, UFS) = prob_parm.rho_r;
    }
}

/**
 * \brief Fill external boundary conditions for ghost cells.
 * 
 * @param x         ghost cell cooridinates.
 * @param s_int     flow state inside of the domain.
 * @param s_ext     flow state to be filled.
 * @param idir      direction (0: x, 1: y, 2: z).
 * @param sgn       high or low boundary (1: low, -1: high).
 * @param time      time.
 * @param geomdata  domain geometry data.
 * @param prob_parm ProbParm data as defined in prob_parm.H and initialised in amrex_probinit.
 * @sa CnsFillExtDir
 * @sa CnsFillExtDir::operator()
 */
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
bcnormal (const amrex::Real x[AMREX_SPACEDIM],
          const amrex::Real s_int[LEN_STATE],
          amrex::Real s_ext[LEN_STATE],
          const int idir,
          const int sgn, 
          const amrex::Real time,
          amrex::GeometryData const& /*geomdata*/,
          ProbParm const& prob_parm)
{
    if (idir == 0 && sgn == 1) { //xlo
        // post shock conditions
        s_ext[URHO] = prob_parm.rho_l;
        s_ext[UMX] = prob_parm.rho_l * prob_parm.u_l;
        s_ext[UMY] = prob_parm.rho_l * prob_parm.v_l;
        s_ext[UMZ] = 0.0;
        s_ext[UEDEN] = prob_parm.rhoe_l + 0.5 * prob_parm.rho_l * (prob_parm.u_l * prob_parm.u_l + prob_parm.v_l * prob_parm.v_l);
        s_ext[UEINT] = prob_parm.rhoe_l;
        s_ext[UTEMP] = prob_parm.T_l;
        s_ext[UFS] = prob_parm.rho_l;
    } else if (idir == 1) {
        if (sgn == 1) { //ylo
            if (x[0] < 1./6.) {
                // FOExtrap
                for (int i = 0; i < LEN_STATE; ++i) s_ext[i] = s_int[i];
            } else {
                // SlipWall
                s_ext[URHO] = s_int[URHO];
                s_ext[UMX] = s_int[UMX];
                s_ext[UMY] =-s_int[UMY];
                s_ext[UMZ] = s_int[UMZ];
                s_ext[UEDEN] = s_int[UEDEN];
                s_ext[UEINT] = s_int[UEINT];
                s_ext[UTEMP] = s_int[UTEMP];
                s_ext[UFS] = s_int[UFS];
            }
        } else { //yhi
            // moving shock conditions
            amrex::Real us = 10. * prob_parm.c_r / sqrt(3.) * 2.; //shock speed
            if (x[0] - us*time < 1./6. + 1./sqrt(3.)*x[1]) {
                s_ext[URHO] = prob_parm.rho_l;
                s_ext[UMX] = prob_parm.rho_l * prob_parm.u_l;
                s_ext[UMY] = prob_parm.rho_l * prob_parm.v_l;
                s_ext[UMZ] = 0.0;
                s_ext[UEDEN] = prob_parm.rhoe_l + 0.5 * prob_parm.rho_l * (prob_parm.u_l * prob_parm.u_l + prob_parm.v_l * prob_parm.v_l);
                s_ext[UEINT] = prob_parm.rhoe_l;
                s_ext[UTEMP] = prob_parm.T_l;
                s_ext[UFS] = prob_parm.rho_l;
            } else {
                s_ext[URHO] = prob_parm.rho_r;
                s_ext[UMX] = 0.0;
                s_ext[UMY] = 0.0;
                s_ext[UMZ] = 0.0;
                s_ext[UEDEN] = prob_parm.rhoe_r;
                s_ext[UEINT] = prob_parm.rhoe_r;
                s_ext[UTEMP] = prob_parm.T_r;
                s_ext[UFS] = prob_parm.rho_r;
            }  
        }
    }
}

/**
 * \brief Modify state data and/or add turbulence to fields after restart.
 * 
 * @param i         x position.
 * @param j         y position.
 * @param k         z position.
 * @param state     state data.
 * @param geomdata  domain geometry data.
 * @param parm      Parm data defined in parm.H.
 * @param prob_parm ProbParm data as defined in prob_parm.H and initialised in amrex_probinit.
 */
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
cns_prob_post_restart (int i, int j, int k, 
                       amrex::Array4<amrex::Real> const& state,
                       amrex::GeometryData const& geomdata,
                       Parm const& parm,
                       ProbParm const& prob_parm)
{
}

#endif