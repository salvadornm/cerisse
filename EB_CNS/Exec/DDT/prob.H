#ifndef CNS_PROB_H
#define CNS_PROB_H

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_ParmParse.H>
#include <AMReX_REAL.H>

#include "CNS.H"
#include "PelePhysics.H"
#include "prob_parm.H"

// Declare your own geometry here and define the build method in prob.cpp
#include "custom_geometry.H"
class GraVent : public CustomGeometry::Register<GraVent>
{
public:
  static const std::string identifier() { return "GraVent"; }

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
              amrex::GeometryData const& geomdata, Parm const& /*parm*/,
              ProbParm const& prob_parm)
{
  // Geometry
  const amrex::Real* prob_lo = geomdata.ProbLo();
  const amrex::Real* dx = geomdata.CellSize();
  const amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];
  const amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1];
  const amrex::Real z = prob_lo[2] + (k + 0.5) * dx[2];
  const amrex::Real rsqr = AMREX_D_TERM(x * x, +y * y, +z * z);

  // Parse mixture details
  amrex::Real volH2;
  {
    amrex::ParmParse pp("prob");
    pp.get("volH2", volH2);
  }

  amrex::GpuArray<amrex::Real, NUM_SPECIES> X = {0.0};
  X[H2_ID] = volH2; // Note: mole frac and vol frac are the same for perfect gases
  X[N2_ID] = (1.0 - volH2) * 0.79;
  X[O2_ID] = (1.0 - volH2) * 0.21;

  // Atmospheric conditions
  amrex::Real p = 1013250.0, Y[NUM_SPECIES] = {0.0};
  amrex::Real T = (rsqr >= 1.0) ? 295.0 : 2500.0; // ignition kernel of 1cm radius

  // Convert to rho, ei, Y
  amrex::Real rho, ei;
  auto eos = pele::physics::PhysicsType::eos();
  eos.X2Y(X.begin(), Y);
  eos.PYT2RE(p, Y, T, rho, ei);

  for (int nf = 0; nf <= NUM_FIELD; ++nf) {
    state(i, j, k, nf * NVAR + URHO) = rho;
    state(i, j, k, nf * NVAR + UMX) = 0.0;
    state(i, j, k, nf * NVAR + UMY) = 0.0;
    state(i, j, k, nf * NVAR + UMZ) = 0.0;
    state(i, j, k, nf * NVAR + UEDEN) = rho * ei;
    for (int n = 0; n < NUM_SPECIES; ++n)
      state(i, j, k, nf * NVAR + UFS + n) = rho * Y[n];
  }
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
bcnormal(const amrex::Real x[AMREX_SPACEDIM], const amrex::Real s_int[LEN_STATE],
         const amrex::Real s_refl[LEN_STATE], amrex::Real s_ext[LEN_STATE],
         const int idir, const int sgn, const amrex::Real time,
         amrex::GeometryData const& /*geomdata*/, ProbParm const& prob_parm)
{
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void prob_post_restart(
  int i, int j, int k, amrex::Array4<amrex::Real> const& state,
  amrex::GeometryData const& geomdata, Parm const& parm, ProbParm const& prob_parm)
{
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void prob_post_timestep(
  int i, int j, int k, const amrex::Real /*curtime*/, const amrex::Real dtlev,
  amrex::Array4<amrex::Real> const& state,
  amrex::Array4<const amrex::Real> const& reactions,
  amrex::GeometryData const& /*geomdata*/, Parm const& /*parm*/, ProbParm const& pp)
{
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void prob_tag_error(
  int i, int j, int k, amrex::Array4<char> const& tagarr,
  amrex::Array4<const amrex::Real> const& sarr, int /*level*/, char tagval,
  const amrex::Real time, amrex::GeometryData const& geomdata, Parm const& /*parm*/,
  ProbParm const& /*pp*/)
{
}

#endif