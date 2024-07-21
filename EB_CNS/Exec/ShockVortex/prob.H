#ifndef CNS_PROB_H
#define CNS_PROB_H

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_REAL.H>
#include <PelePhysics.H>

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
 * @param prob_parm ProbParm data as defined in prob_parm.H and initialised in
 * amrex_probinit.
 */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
prob_initdata(int i, int j, int k, amrex::Array4<amrex::Real> const& state,
              amrex::GeometryData const& geomdata, ProbParm const& prob_parm)
{
  // Geometry
  const amrex::Real* prob_lo = geomdata.ProbLo();
  const amrex::Real* dx = geomdata.CellSize();
  const amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];
  const amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1];

  const amrex::Real R =
    pele::physics::Constants::RU / pele::physics::Constants::AIRMW;
  const amrex::Real g = pele::physics::Constants::gamma;

  // Vortex location
  const amrex::Real x_c = 0.25;
  const amrex::Real y_c = 0.5;
  // Vortex size
  const amrex::Real a = 0.075;
  const amrex::Real b = 0.175;
  // Vortex strength
  const amrex::Real M_v = 0.9;
  const amrex::Real v_m = M_v * sqrt(g);
  // Distance from vortex
  const amrex::Real r = sqrt((x - x_c) * (x - x_c) + (y - y_c) * (y - y_c));

  // Find states at the position
  amrex::Real rho = x <= 0.5 ? prob_parm.rho_l : prob_parm.rho_r;
  amrex::Real u = x <= 0.5 ? prob_parm.u_l : prob_parm.u_r;
  amrex::Real v = 0.0;
  amrex::Real rho_e = (x <= 0.5 ? prob_parm.p_l : prob_parm.p_r) / (g - 1.0);

  if (r <= b) {
    const amrex::Real sin_theta = (y - y_c) / r;
    const amrex::Real cos_theta = (x - x_c) / r;

    amrex::Real t_u = prob_parm.p_l / (prob_parm.rho_l * R); // vortex is upstream
    amrex::Real radial_term, temperature;

    if (r <= a) {
      u -= v_m * r / a * sin_theta;
      v += v_m * r / a * cos_theta;

      radial_term = -2.0 * b * b * log(b) - (0.5 * a * a) + (2.0 * b * b * log(a)) +
                    (0.5 * b * b * b * b / (a * a));
      amrex::Real t_a =
        t_u - (g - 1.0) * pow(v_m * a / (a * a - b * b), 2) * radial_term / (R * g);
      radial_term = 0.5 * (1.0 - r * r / (a * a));
      temperature = t_a - (g - 1.0) * v_m * v_m * radial_term / (R * g);
    } else {
      u -= v_m * a * (r - b * b / r) / (a * a - b * b) * sin_theta;
      v += v_m * a * (r - b * b / r) / (a * a - b * b) * cos_theta;

      radial_term = -2.0 * b * b * log(b) - (0.5 * r * r) + (2.0 * b * b * log(r)) +
                    (0.5 * b * b * b * b / (r * r));
      temperature =
        t_u - (g - 1.0) * pow(v_m * a / (a * a - b * b), 2) * radial_term / (R * g);
    }
    rho_e = prob_parm.p_l * pow(temperature / t_u, g / (g - 1.0)) / (g - 1.0);
  }

  // Finally, set the states
  for (int nf = 0; nf <= NUM_FIELD; ++nf) {
    state(i, j, k, nf * NVAR + URHO) = rho;
    state(i, j, k, nf * NVAR + UMX) = rho * u;
    state(i, j, k, nf * NVAR + UMY) = rho * v;
    state(i, j, k, nf * NVAR + UMZ) = 0.0;
    state(i, j, k, nf * NVAR + UEDEN) = rho_e + 0.5 * rho * (u * u + v * v);
    state(i, j, k, nf * NVAR + UFS) = rho;
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
    for (int nf = 0; nf <= NUM_FIELD; ++nf) {
      s_ext[nf * NVAR + URHO] = prob_parm.rho_l;
      s_ext[nf * NVAR + UMX] = prob_parm.rho_l * prob_parm.u_l;
      s_ext[nf * NVAR + UMY] = 0.0;
      s_ext[nf * NVAR + UMZ] = 0.0;
      s_ext[nf * NVAR + UEDEN] =
        prob_parm.rhoe_l + 0.5 * prob_parm.rho_l * prob_parm.u_l * prob_parm.u_l;
      s_ext[nf * NVAR + UFS] = prob_parm.rho_l;
    }
  }
}

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