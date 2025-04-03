#ifndef CNS_PROB_H
#define CNS_PROB_H

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_REAL.H>

#include "CNS.H"
#include "index_macros.H"
#include "mechanism.H"
#include "prob_parm.H"

#include "custom_geometry.H"  

class Combustor_NTNU : public CustomGeometry::Register<Combustor_NTNU>
{
public:
  static const std::string identifier() { return "Combustor_NTNU"; }

  void build(const Geometry& geom, const int max_coarsening_level) override;
};

/**
 * \brief Initialise state data.
 *
 * @param i,j,k    x, y, z index.
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
  auto eos = pele::physics::PhysicsType::eos();
  amrex::Real rho, eint;
  eos.PYT2RE(prob_parm.interior_p, (amrex::Real*)prob_parm.massfrac.begin(),
             prob_parm.interior_T, rho, eint);

  state(i, j, k, URHO) = rho;
  state(i, j, k, UMX) = amrex::Real(0.0);
  state(i, j, k, UMY) = amrex::Real(0.0);
  state(i, j, k, UMZ) = amrex::Real(0.0);
  state(i, j, k, UEDEN) = rho * eint; 
  for (int n = 0; n < NUM_SPECIES; ++n) {
    state(i, j, k, UFS + n) = state(i, j, k, URHO) * prob_parm.massfrac[n];
  }
}

/**
 * \brief Fill external boundary conditions for ghost cells.
 *
 * @param x         ghost cell cooridinates.
 * @param s_int     flow state inside of the domain.
 * @param s_ext     flow state to be filled.
 * @param idir      direction, 0:x, 1:y, 2:z.
 * @param sgn       high or low boundary (1: low, -1: high).
 * @param time      time.
 * @param geomdata  domain geometry data.
 * @param prob_parm ProbParm data as defined in prob_parm.H and initialised in
 * amrex_probinit.
 * @sa CnsFillExtDir
 * @sa CnsFillExtDir::operator()
 */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void bcnormal(
  const amrex::Real x[AMREX_SPACEDIM], const amrex::Real* /*s_int[LEN_STATE]*/,
  const amrex::Real* /*s_refl[LEN_STATE]*/, amrex::Real s_ext[LEN_STATE],
  const int idir, const int sgn, const amrex::Real /*time*/,
  amrex::GeometryData const& /*geomdata*/, ProbParm const& prob_parm)
{
  // Z=0  low boundary 
  if (idir == 2 && sgn == 1) {
    for (int nf = 0; nf <= NUM_FIELD; ++nf) {
      for (int nc = 0; nc < NVAR; ++nc) {
        s_ext[nf * NVAR + nc] = prob_parm.inflow_state[nc];
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

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
prob_get_aux_name(amrex::Vector<std::string>& aux_name)
{
}

#endif