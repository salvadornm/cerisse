#ifndef PIPE_IBM_PROB_H_
#define PIPE_IBM_PROB_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_ParmParse.H>
#include <PelePhysics.H>
#include <ReactorBase.H>
#include "Closures.h"
#include "RHS.h"
#include "LES.h"
#include <ibm_solver.h>

using namespace amrex;

namespace PROB {

inline constexpr int ibm_eorder = 1;

// Match the WALE LES closure used by the CamHc1 IBM case.
struct LESparm {
  static constexpr Real Cs = 0.17;
  static constexpr int order = 2;
  static constexpr Real Scsgs = 0.4;
  static constexpr Real Pr_o_Prsgs = 1.0;
  static constexpr bool fixDelta = false;
};

using ProbClosures = closures_dt<indicies_t, transport_Pele_t,
    multispecies_pele_gas_t<indicies_t>, WALE_t<LESparm,indicies_t>>;

struct ProbParm {
  static constexpr Real P0 = 101325.0;
  static constexpr Real T0 = 300.0;
  static constexpr Real Uin = 70.0;
  static constexpr Real pipe_radius = 0.05;
  static constexpr Real z_in = 0.0;
  static constexpr Real z_out = 2.5;
  static constexpr Real port_tol = 1.0e-5;

  Real rho = 0.0;
  Real eint = 0.0;
  Real Y[NUM_SPECIES] = {};

  ProbParm() {
    Y[O2_ID] = Real(0.233);
    Y[N2_ID] = Real(0.767);
    ProbClosures cls;
    cls.PYT2R(P0, Y, T0, rho);
    cls.RYP2E(rho, Y, P0, eint);
  }
};

struct skewparm_t {
  static constexpr bool dissipation = true;
  static constexpr int order = 4;
  static constexpr Real C2skew = 0.5;
  static constexpr Real C4skew = 0.016;
};

struct viscous_param_t {
  static constexpr bool ATF = false;
  static constexpr Real thickfactor = 1.0;
  static constexpr int order = 2;
  static constexpr bool use_LES = true;
};

struct ibmparm_t {
  static constexpr int interp_order = 1;
  static constexpr int extrap_order = ibm_eorder;
  static constexpr Real alpha = 1.5;
  static constexpr int interp_order_surf = 1;
  static constexpr int extrap_order_surf = 1;
  static constexpr Real alpha_surf = 0.6;
  static constexpr int ghost_layers = 1;
  // The watertight cylinder bounds the pipe fluid; its exterior is solid.
  static constexpr bool interior_is_solid = false;
};

using ProbRHS = rhs_dt<skew_t<skewparm_t,ProbClosures>,
    viscousLES_t<viscous_param_t,ProbClosures>, no_source_t>;

template <typename param, typename cls_t>
class pipe_ibm_wall_t {
public:
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  static void compute_surfIB(
      const Array1D<Real,0,AMREX_SPACEDIM-1>& xyz,
      const Array1D<Real,0,AMREX_SPACEDIM-1>& norm,
      const Array1D<Real,0,AMREX_SPACEDIM-1>& t1,
      const Array1D<Real,0,AMREX_SPACEDIM-1>& t2,
      Array2D<Real,0,ibm_eorder+1,0,cls_t::NPRIM-1>& q,
      int /*type_solid_bc*/, const cls_t* /*cls*/) {
    const Real r2 = xyz(0)*xyz(0) + xyz(1)*xyz(1);
    const bool cap_normal = amrex::Math::abs(norm(2)) > Real(0.9);
    const bool inlet = cap_normal &&
        amrex::Math::abs(xyz(2)-param::z_in) < param::port_tol &&
        r2 <= param::pipe_radius*param::pipe_radius;
    const bool outlet = cap_normal &&
        amrex::Math::abs(xyz(2)-param::z_out) < param::port_tol &&
        r2 <= param::pipe_radius*param::pipe_radius;

    // Default pipe wall: adiabatic, no-slip, impermeable, zero-gradient p.
    q(1,cls_t::QU) = Real(0.0);
    q(1,cls_t::QV) = Real(0.0);
    q(1,cls_t::QW) = Real(0.0);
    q(1,cls_t::QPRES) = q(2,cls_t::QPRES);
    q(1,cls_t::QT) = q(2,cls_t::QT);
    for (int sp = 0; sp < NUM_SPECIES; ++sp)
      q(1,cls_t::QFS+sp) = q(2,cls_t::QFS+sp);

    if (inlet) {
      // interior_is_solid=false flips the bottom-cap normal into +z.
      q(1,cls_t::QU) = param::Uin;
      q(1,cls_t::QPRES) = param::P0;
      q(1,cls_t::QT) = param::T0;
      for (int sp = 0; sp < NUM_SPECIES; ++sp)
        q(1,cls_t::QFS+sp) = Real(0.0);
      q(1,cls_t::QFS+O2_ID) = Real(0.233);
      q(1,cls_t::QFS+N2_ID) = Real(0.767);
    } else if (outlet) {
      q(1,cls_t::QU) = q(2,cls_t::QU);
      q(1,cls_t::QV) = q(2,cls_t::QV);
      q(1,cls_t::QW) = q(2,cls_t::QW);
      q(1,cls_t::QPRES) = param::P0;
      q(1,cls_t::QT) = q(2,cls_t::QT);
    }
    amrex::ignore_unused(t1,t2);
  }
};

using TypeWall = pipe_ibm_wall_t<ProbParm,ProbClosures>;
using ProbIB = ibm_solver_t<TypeWall,ibmparm_t,ProbClosures>;

inline void inputs() {}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void prob_initdata(int i, int j, int k, Array4<Real> const& state,
    GeometryData const&, ProbClosures const& cls, ProbParm const& pp) {
  state(i,j,k,cls.URHO) = pp.rho;
  state(i,j,k,cls.UMX) = 0.0;
  state(i,j,k,cls.UMY) = 0.0;
  state(i,j,k,cls.UMZ) = pp.rho*ProbParm::Uin;
  state(i,j,k,cls.UET) = pp.rho*pp.eint +
      Real(0.5)*pp.rho*ProbParm::Uin*ProbParm::Uin;
  for (int sp = 0; sp < NUM_SPECIES; ++sp)
    state(i,j,k,cls.UFS+sp) = pp.rho*pp.Y[sp];
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void bcnormal(const Real[AMREX_SPACEDIM], Real,
    const Real s_int[ProbClosures::NCONS], const Real[ProbClosures::NCONS],
    Real s_ext[ProbClosures::NCONS], int, int, Real, GeometryData const&,
    ProbClosures const&, ProbParm const&) {
  for (int n = 0; n < ProbClosures::NCONS; ++n) s_ext[n] = s_int[n];
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void user_tagging(int i, int j, int k, int, auto& tags, const auto&,
    const auto&, const auto& geom, const ProbParm&, int level) {
  const auto* lo = geom.ProbLo();
  const auto* dx = geom.CellSize();
  const Real x = lo[0] + (i+Real(0.5))*dx[0];
  const Real y = lo[1] + (j+Real(0.5))*dx[1];
  const Real z = lo[2] + (k+Real(0.5))*dx[2];
  // Refine the pipe interior and a thin band around its wall at level 1.
  tags(i,j,k) = level == 0 && x*x + y*y < Real(0.055)*Real(0.055) &&
      z > Real(-0.02) && z < Real(2.52);
}

} // namespace PROB
#endif
