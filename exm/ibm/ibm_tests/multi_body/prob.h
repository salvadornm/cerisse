#ifndef CNS_PROB_H_
#define CNS_PROB_H_

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_ParmParse.H>
#include <AMReX_AmrLevel.H>
#include <Closures.h>
#include <RHS.h>
#include <ibm_solver.h>
#include <Constants.h>
#include <ibm_wallmodel.h>

using namespace amrex;
using namespace universal_constants;

namespace PROB {

// Physical constants
static constexpr Real Mach     = 4.0;
static constexpr Real Mw       = 28.96e-3;
static constexpr Real gam      = 1.4;
static constexpr Real Rgas     = gas_constant / Mw;
static constexpr Real Reynolds = 10000;
static constexpr Real Pr       = 0.7;
static constexpr Real Cv       = Rgas / (gam - 1.0);
static constexpr Real Cp       = gam * Cv;
static constexpr Real viscos   = 1.0 / Reynolds;
static constexpr Real lambda   = viscos * Cp / Pr;

//////////////////////////// Physical modelling ////////////////////////////////
struct ProbParm
{
  // Free-stream conditions
  Real p_oo    = 100000.0;
  Real T_oo    = 300.0;
  Real rho_oo  = p_oo / (Rgas * T_oo);
  Real c_oo    = sqrt(gam * Rgas * T_oo);
  Real u_oo    = c_oo * Mach;
  Real eint_oo = rho_oo * Cv * T_oo;
  Real kin_oo  = 0.5 * rho_oo * u_oo * u_oo;

  // Right state (quiescent)
  Real p_r     = p_oo;
  Real T_r     = T_oo;
  Real rho_r   = p_r / (Rgas * T_r);
  Real u_r     = 0.0;
  Real eint_r  = rho_r * Cv * T_r;

  // Initial shock position
  Real xshock = -1.5;
};

// Viscous solver parameters
struct methodparm_t {
  static constexpr int  order = 2;
  static constexpr Real conductivity = lambda;
  static constexpr Real viscosity    = viscos;
  static constexpr bool use_LES = false;
};

// Skew-symmetric method parameters (unused but kept for compatibility)
struct skewparm_t {
  static constexpr bool dissipation = true;
  static constexpr int  order = 4;
  static constexpr Real C2skew = 1.5, C4skew = 0.016;
};

// IBM parameters — isothermal wall for heat flux computation
struct ibmparm_t {
  static constexpr int  interp_order = 1;
  static constexpr int  extrap_order = 1;
  static constexpr Real alpha = 0.6;

  static constexpr int  interp_order_surf = 1;
  static constexpr int  extrap_order_surf = 1;
  static constexpr Real alpha_surf = 0.6;

  static constexpr int  ghost_layers = 1;
  static constexpr bool interior_is_solid = true;

  // Wall temperature for isothermal BC (enables nonzero dTdn / heat flux)
  static constexpr Real Twall = 300.0;
};

// CLOSURES
typedef closures_dt<indicies_t, transport_const_t<methodparm_t>,
                    calorifically_perfect_gas_t<indicies_t>> ProbClosures;

// NUMERICAL SCHEME: Riemann + viscous (needed for shear stress and heat flux)
typedef rhs_dt<riemann_t<false, ProbClosures>, viscous_t<methodparm_t, ProbClosures>, no_source_t> ProbRHS;

// IBM wall type: isothermal no-slip → nonzero dTdn (heat flux) and tau (shear)
typedef ibm_isothermal_noslip_wall_t<ibmparm_t, ProbClosures> TypeWall;
typedef ibm_solver_t<TypeWall, ibmparm_t, ProbClosures> ProbIB;

// Static geometry: no update needed
inline void update_geometry(Real /*time*/,
                            Vector<GeomType>& /*geom_a*/,
                            int /*ngeom*/) {}

void inline inputs() {
  amrex::Print() << " ****** Starting ... ******* " << std::endl;
  amrex::Print() << " Mach 4 flow over 7 static cylinders (IBM)" << std::endl;
}

//////////////////////////// Initial conditions ////////////////////////////////
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void prob_initdata(int i, int j, int k, amrex::Array4<amrex::Real> const& state,
    amrex::GeometryData const& geomdata, ProbClosures const& cls, ProbParm const& pparm) {

  const Real* prob_lo = geomdata.ProbLo();
  const Real* dx      = geomdata.CellSize();
  Real x = prob_lo[0] + (i + 0.5_rt) * dx[0];

  Real rhot, eint, u[3] = {0.0};

  if (x < pparm.xshock) {
    rhot = pparm.rho_oo;
    u[0] = pparm.u_oo;
    eint = pparm.eint_oo;
  } else {
    rhot = pparm.rho_r;
    u[0] = pparm.u_r;
    eint = pparm.eint_r;
  }

  state(i, j, k, cls.URHO) = rhot;
  state(i, j, k, cls.UMX)  = rhot * u[0];
  state(i, j, k, cls.UMY)  = Real(0.0);
  state(i, j, k, cls.UMZ)  = Real(0.0);
  state(i, j, k, cls.UET)  = eint + Real(0.5) * rhot * u[0] * u[0];
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void user_tagging(int i, int j, int k, int /*nt*/, auto& tagfab, const auto& sdatafab,
                  const Array4<uint8_t>& ibfab, const auto& /*geomdata*/,
                  const ProbParm& /*pparm*/, int level) {
  // Density-gradient based tagging
  int URHO = ProbClosures::URHO;
  Real rhop  = sdatafab(i, j, k, URHO);
  Real drhox = std::abs(sdatafab(i+1, j, k, URHO) - sdatafab(i-1, j, k, URHO));
  Real drhoy = std::abs(sdatafab(i, j+1, k, URHO) - sdatafab(i, j-1, k, URHO));
  Real rhofluc = std::sqrt(drhox * drhox + drhoy * drhoy) / rhop;
  Real threshold[3] = {0.3_rt, 0.6_rt, 1000_rt};
  tagfab(i, j, k) = (rhofluc > threshold[level]);

  // Refine near IBM body: tag 3x3 stencil around every ghost cell
  if (ibfab(i, j, k, 1)) {
    for (int ii = -1; ii <= 1; ii++) {
      for (int jj = -1; jj <= 1; jj++) {
        tagfab(i + ii, j + jj, k) = true;
      }
    }
  }
}

//////////////////////////// Boundary conditions ///////////////////////////////
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
bcnormal(const amrex::Real x[AMREX_SPACEDIM], amrex::Real dratio,
         const amrex::Real s_int[ProbClosures::NCONS],
         const amrex::Real s_refl[ProbClosures::NCONS],
         amrex::Real s_ext[ProbClosures::NCONS],
         const int idir, const int sgn, const amrex::Real time,
         amrex::GeometryData const& /*geomdata*/,
         ProbClosures const& closures, ProbParm const& pparm)
{
  const int URHO = ProbClosures::URHO;
  const int UMX  = ProbClosures::UMX;
  const int UMY  = ProbClosures::UMY;
  const int UMZ  = ProbClosures::UMZ;
  const int UET  = ProbClosures::UET;
  const int face = (idir + 1) * sgn;

  switch (face)
  {
    case 1:  // x-lo — inflow
      s_ext[URHO] = pparm.rho_oo;
      s_ext[UMX]  = pparm.rho_oo * pparm.u_oo;
      s_ext[UMY]  = 0.0;
      s_ext[UMZ]  = 0.0;
      s_ext[UET]  = pparm.eint_oo + pparm.kin_oo;
      break;
    default:  // outflow: zero-gradient (s_ext already filled with s_int)
      break;
  }
}

} // namespace PROB
#endif
