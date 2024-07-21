#ifndef CNS_PROB_H
#define CNS_PROB_H

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_REAL.H>

#include "CNS.H"
#include "custom_geometry.H"
#include "prob_parm.H"

// Declare your own geometry here and define the build method in prob.cpp
class DoubleRamp : public CustomGeometry::Register<DoubleRamp>
{
public:
  static const std::string identifier() { return "double_ramp"; }

  void build(const Geometry& geom, const int max_coarsening_level) override;
};

/**
 * \brief Initialise state data.
 *
 * @param i         x position.
 * @param j         y position.
 * @param k         z position.
 * @param state     state data.
 * @param geomdata  domain geometry data.
 * @param parm      Parm data defined in parm.H.
 * @param prob_parm ProbParm data as defined in prob_parm.H and initialised in
 * amrex_probinit.
 */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
prob_initdata(int i, int j, int k, amrex::Array4<amrex::Real> const& state,
              amrex::GeometryData const& geomdata, ProbParm const& prob_parm)
{
  // Geometry
  // const amrex::Real* prob_lo = geomdata.ProbLo();
  // const amrex::Real* dx = geomdata.CellSize();
  // const amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];
  // const amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1];

  // Set the states
  // if ((x < 1.0 /*&& y > 0.2/**/) /*|| (x <= 0.5-y && y <= 0.2)*/) {
  state(i, j, k, URHO) = prob_parm.rho;
  state(i, j, k, UMX) = prob_parm.rho * prob_parm.u;
  state(i, j, k, UMY) = 0.0;
  state(i, j, k, UMZ) = 0.0;
  state(i, j, k, UEDEN) =
    prob_parm.rhoe + 0.5 * prob_parm.rho * prob_parm.u * prob_parm.u;
  for (int n = 0; n < NUM_SPECIES; ++n)
    state(i, j, k, UFS + n) = prob_parm.rho * prob_parm.massfrac[n];
  /*} else if (x > 0.5-y && y <= 0.2) {
      state(i, j, k, URHO) = prob_parm.rho_l;
      state(i, j, k, UMX) = 0.0;
      state(i, j, k, UMY) = 0.25 * prob_parm.rho_l * prob_parm.u_l;
      state(i, j, k, UMZ) = 0.0;
      state(i, j, k, UEDEN) = prob_parm.rhoe_l + 0.5 * 0.125 * prob_parm.rho_l *
     prob_parm.u_l * prob_parm.u_l; state(i, j, k, UFS) = prob_parm.rho_l;*/
  // } else {
  //     state(i, j, k, URHO) = prob_parm.rho_r;
  //     state(i, j, k, UMX) = 0.0;
  //     state(i, j, k, UMY) = 0.0;
  //     state(i, j, k, UMZ) = 0.0;
  //     state(i, j, k, UEDEN) = prob_parm.rhoe_r;
  //     state(i, j, k, UFS) = prob_parm.rho_r;
  // }
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
 * @param prob_parm ProbParm data as defined in prob_parm.H and initialised in
 * amrex_probinit.
 * @sa CnsFillExtDir
 * @sa CnsFillExtDir::operator()
 */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
bcnormal(const amrex::Real x[AMREX_SPACEDIM], const amrex::Real s_int[LEN_STATE],
         const amrex::Real s_refl[LEN_STATE], amrex::Real s_ext[LEN_STATE],
         const int idir, const int sgn, const amrex::Real time,
         amrex::GeometryData const& /*geomdata*/, ProbParm const& prob_parm)
{
  if (idir == 0 && sgn == 1) { // xlo
    // if (x[1] > 0.5) {
    // inlet: post shock conditions
    s_ext[URHO] = prob_parm.rho;
    s_ext[UMX] = prob_parm.rho * prob_parm.u;
    s_ext[UMY] = 0.0;
    s_ext[UMZ] = 0.0;
    s_ext[UEDEN] = prob_parm.rhoe + 0.5 * prob_parm.rho * prob_parm.u * prob_parm.u;
    for (int n = 0; n < NUM_SPECIES; ++n)
      s_ext[UFS + n] = prob_parm.rho * prob_parm.massfrac[n];
    // } else {
    //     s_ext[URHO] = prob_parm.rho_r;
    //     s_ext[UMX] = 0.0;
    //     s_ext[UMY] = 0.0;
    //     s_ext[UMZ] = 0.0;
    //     s_ext[UEDEN] = prob_parm.rhoe_r;
    //     s_ext[UFS] = prob_parm.rho_r;
    // }
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
 * @param prob_parm ProbParm data as defined in prob_parm.H and initialised in
 * amrex_probinit.
 */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
prob_post_restart(int i, int j, int k, amrex::Array4<amrex::Real> const& state,
                  amrex::GeometryData const& geomdata, ProbParm const& prob_parm)
{
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
prob_post_timestep(int i, int j, int k, const amrex::Real /*curtime*/,
                   const amrex::Real dtlev, amrex::Array4<amrex::Real> const& state,
                   amrex::Array4<const amrex::Real> const& reactions,
                   amrex::GeometryData const& /*geomdata*/, ProbParm const& pp)
{
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
prob_tag_error(int i, int j, int k, amrex::Array4<char> const& tagarr,
               amrex::Array4<const amrex::Real> const& /*sarr*/, int level,
               char tagval, const amrex::Real /*time*/,
               amrex::GeometryData const& geomdata, ProbParm const& /*pp*/)
{
}

#endif