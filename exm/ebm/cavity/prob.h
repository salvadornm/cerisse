// Two-dimensional laminar flow over an open rectangular cavity.
//
// This setup targets case 2M6 (also called L2) of Rowley, Colonius &
// Basu, JFM 455 (2002): M_inf = 0.6, L/D = 2, L/theta_0 = 52.8 and
// Re_theta = 56.8.  Lengths and the embedded-boundary geometry are defined
// in SI units; see eb_geometry.cpp for D and L.

#ifndef CNS_PROB_H_
#define CNS_PROB_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_ParmParse.H>

#include "Closures.h"
#include "RHS.h"

#if CNS_USE_EB    
#include <ebm.h>
#include <walltypes.h>
#endif

#include <bc_types.h>
#include <cmath>

using namespace amrex;

namespace PROB {

struct CavityParm {
  // Free-stream air at one atmosphere.
  static constexpr Real rho_inf = 1.177;          // kg/m^3
  static constexpr Real T_inf = 300.0;            // K
  static constexpr Real viscosity = 1.845e-5;     // Pa s
  static constexpr Real conductivity = 0.0262;    // W/(m K)
  static constexpr Real R = 287.05;               // J/(kg K)
  static constexpr Real gamma = 1.4;

  // Paper parameters for run 2M6/L2.
  static constexpr Real Ma = 0.6;
  static constexpr Real Re_theta = 56.8;
  static constexpr Real theta_0 = 4.274e-6;        // m, at cavity leading edge

  // Virtual origin used by the polynomial boundary-layer approximation.
  // With the leading edge at x = 3.9 D, this gives theta_0 ~= 4.28 um.
  static constexpr Real x_virtual = 1.900e-4;      // m
};

// LES closures
struct LESparm {
  static constexpr Real Cs = 0.0;           // zero disables the SGS model
  static constexpr int order = 2;
  static constexpr Real Scsgs = 0.4;
  static constexpr Real Pr_o_Prsgs = 0.9;
  static constexpr bool fixDelta = false;
};


typedef closures_dt<indicies_stat_t, transport_const_t<CavityParm>, calorifically_perfect_gas_t<indicies_t>, Smagorinsky_t<LESparm,indicies_t>>ProbClosures;

// Dimensional reference state shared by initialization and boundary conditions.
struct ProbParm {
  ProbClosures pp_pc;

  const Real T_0 = CavityParm::T_inf;
  const Real p_0 = pres_atm2si;
  Real rho_0, eint_0, c_0;
  Real vel_inlet[3] = {0.0, 0.0, 0.0};
  Real mass_flux;

  ProbParm() {
    pp_pc.PYT2R(p_0, nullptr, T_0, rho_0);
    pp_pc.RYP2E(rho_0, nullptr, p_0, eint_0);
    pp_pc.PYT2Cs(p_0, nullptr, T_0, c_0);
    vel_inlet[0] = CavityParm::Ma * c_0;

    // Re_theta = rho_inf U_inf theta_0 / mu, so rho_inf U_inf is
    // prescribed as a mass flux at the inflow boundary.
    mass_flux = CavityParm::Re_theta * CavityParm::viscosity / CavityParm::theta_0;
  }
};

// numerical method parameters 
struct skewparm_t {

  public:

  static constexpr bool dissipation = false; // paper relies on physical viscosity
  static constexpr int order = 4;            // supported values: 2, 4 or 6
  static constexpr Real C2skew = 0.1;
  static constexpr Real C4skew = 0.016;
};


struct viscous_param_t {

  public:
  static constexpr int order = 4;
  static constexpr bool use_LES = false; // the paper is a two-dimensional DNS
};

struct wall_param {

  public:

  // The paper uses a no-slip, isothermal wall at the free-stream temperature.
  static constexpr Real Twall = CavityParm::T_inf;
  static constexpr bool solve_diffwall = true;

};


template <typename cls_t > class user_source_t;

// Sixth-order skew-symmetric convection, molecular diffusion and a pressure
// sponge near the three far-field boundaries.  No LES or chemistry is used.
typedef rhs_dt<skew_t<skewparm_t, ProbClosures>, viscous_t<viscous_param_t, ProbClosures>, user_source_t<ProbClosures>> ProbRHS;

//typedef rhs_dt<skew_t<skewparm_t, ProbClosures>, viscousLES_t<user_source_t<ProbClosures>, ProbClosures>, user_source_t<ProbClosures>> ProbRHS;


#if CNS_USE_EB    
typedef isothermal_wall_t<wall_param,ProbClosures> TypeWall;
typedef ebm_t<TypeWall,wall_param,ProbClosures> ProbEB;
#endif

typedef manual_bc_t<ProbClosures> GlobalBC;


void inline inputs() {
  //	
}

// initial condition
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
prob_initdata(int i, int j, int k, Array4<Real> const &state,
              GeometryData const &geomdata, ProbClosures const &cls,
              ProbParm const &prob_parm) {
  const Real *prob_lo = geomdata.ProbLo();
  const Real *dx = geomdata.CellSize();

  const Real x = prob_lo[0] + (i + Real(0.5)) * dx[0];
  const Real y = prob_lo[1] + (j + Real(0.5)) * dx[1];

  const Real rho = prob_parm.rho_0;
  const Real eint = prob_parm.eint_0;

  // WARNING: initialize all three stored momentum components even in a 2-D build
  Real u[3] = {prob_parm.vel_inlet[0], 0.0, 0.0};

  // Pohlhausen approximation to the paper's initial Blasius layer:
  //   u/U_inf = 2 eta - 2 eta^3 + eta^4,  eta = y/delta.
  // The profile spans the cavity mouth; cells below y=0 start at rest.
  const Real streamwise_distance = x + CavityParm::x_virtual;
  if (y <= 0.0) {
    u[0] = 0.0;
  } else if (streamwise_distance > 0.0) {
    const Real delta = 5.29 *
      sqrt(CavityParm::viscosity * streamwise_distance /
           (CavityParm::rho_inf * prob_parm.vel_inlet[0]));
    if (y < delta) {
      const Real eta = y / delta;
      u[0] *= 2.0 * eta - 2.0 * eta * eta * eta
              + eta * eta * eta * eta;
    }
  }

  const Real kinetic_energy = Real(0.5) * rho *(u[0] * u[0] + u[1] * u[1] + u[2] * u[2]);
  state(i, j, k, cls.URHO) = rho;
  state(i, j, k, cls.UMX)  = rho * u[0];
  state(i, j, k, cls.UMY)  = rho * u[1];
  state(i, j, k, cls.UMZ)  = rho * u[2];
  state(i, j, k, cls.UET)  = rho * eint + kinetic_energy;
  
}
/////////////////////////////// BC /////////////////////////////////////////////
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
bcnormal(const Real /*x*/[AMREX_SPACEDIM], Real /*dratio*/,
         const Real s_int[ProbClosures::NCONS],
         const Real /*s_refl*/[ProbClosures::NCONS],
         Real s_ext[ProbClosures::NCONS], const int idir,
         const int sgn, const Real /*time*/, GeometryData const & /*geomdata*/,
         ProbClosures const &closures, ProbParm const &prob_parm) {

  const int face = (idir+1)*sgn; // +/-1 (1D) +/- 2 (2D) +/- 3 (3D)

  switch(face)
  {
    case  1:  // west x
	    {                  
      GlobalBC::bc_inlet_fixmassflow(1.0,0.0,0.0,&closures,
        prob_parm.mass_flux, prob_parm.T_0, nullptr, s_int, s_ext);
      break;
      }
    case -1:  // EAST         x  
      GlobalBC::bc_fixP(-1.0,0.0,0.0,&closures,prob_parm.p_0, s_int, s_ext);  
      break;
    case  2:  // SOUTH        y  
      GlobalBC::bc_fixP(0.0,1.0,0.0,&closures,prob_parm.p_0, s_int, s_ext); 
      break;
    case -2:  // NORTH        y 
      GlobalBC::bc_fixP(0.0,-1.0,0.0,&closures,prob_parm.p_0, s_int, s_ext); 
      break;
    default:

      break; 
  }
}
///////////////////////////////AMR//////////////////////////////////////////////
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
user_tagging(int i, int j, int k, int /*nt_level*/, auto &tagfab,
             const auto & /*sdatafab*/, const auto& /*ebflag*/,
             const auto &geomdata, const ProbParm & /*prob_parm*/, int level) {

  const Real *prob_lo = geomdata.ProbLo();
  const Real *dx = geomdata.CellSize();
  Real x = prob_lo[0] + (i + Real(0.5)) * dx[0];
  Real y = prob_lo[1] + (j + Real(0.5)) * dx[1];
  bool refine = false;
        
  switch (level)

    {
      case 0:
        // Level 1: refine the cavity and its near acoustic field.
        //refine = (x > 0.0) && (x < 1.2e-3) && (y < 1.1e-3);
         
        refine = (y < 1.1e-3);

        break;
      case 1:
        refine = (x > 0) && (y < 1.038e-3) && (x < 0.0012972);
        break;
      case 2:
        refine= (x > 4.350e-4) && (x < 6.690e-4) && (y < 0);
        break;
      case 3:
        // refine= (z > 0.035) && (z < 0.07);    
        break;  
        
      default:

      break;
    }

  // refine = true; // temp

  tagfab(i,j,k) = refine;
}
///////////////////////////////SOURCE TERM /////////////////////////////////////
template <typename cls_t>
class user_source_t {
  public:

  // Interface flags required by rhs_dt; chemistry and ATF are disabled.
  bool static constexpr ATF = false; // use adaptive thickening factor
  static constexpr Real thickfactor = 0.0; // thickening factor

  bool static constexpr do_reactions = false;

  static constexpr int order = 2;
  static constexpr bool use_LES = false;
  bool static constexpr mask_cells_boundary = true;

  inline static int src_dt = 0;
  static void inline increase_src_dt(){src_dt++;}
  // Relax pressure in sponge strips adjacent to the inflow, outflow and top
  // boundary. Momentum and energy are changed with density so velocity and
  // specific total energy remain unchanged by the sponge.
  void inline src(const Geometry& geomdata, const amrex::MFIter &mfi,
                  const amrex::Array4<const amrex::Real> &prims,
                  const amrex::Array4<amrex::Real> &rhs, const cls_t *cls_d,
                  amrex::Real /*dt*/, amrex::Real /*real_time*/,
                  const auto& /*ebflag*/){

    const Box& bxg = mfi.tilebox();
    // const Box& bxg = mfi.growntilebox(cls_t::NGHOST);
    const Real *prob_lo = geomdata.ProbLo();
    const Real *dx = geomdata.CellSize();

    ProbParm const prob_parm;
    const auto& cls = *cls_d;

    // Outlet and upper-boundary sponges. The inflow is controlled only by
    // NSCBC, with no volumetric pressure sponge at xlo.
    constexpr Real x0_hi = 1.1145e-3;
    constexpr Real xw_hi = 1.2972e-3;
    constexpr Real y0_hi = 0.998e-3;
    constexpr Real yw_hi = 1.1807e-3;
    constexpr Real inverse_relaxation_time = 2.0e5; // 1/(5 microseconds)

    amrex::ParallelFor(bxg,
      [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
      
      const Real x = prob_lo[0] + (i + Real(0.5)) * dx[0];
      const Real y = prob_lo[1] + (j + Real(0.5)) * dx[1];
      const bool buffer = (x >= x0_hi) || (y >= y0_hi);

      if (buffer) {
        const Real d_xhi = (x >= x0_hi) ? (x - x0_hi) / (xw_hi - x0_hi) : 0.0;
        const Real d_yhi = (y >= y0_hi) ? (y - y0_hi) / (yw_hi - y0_hi) : 0.0;
        const Real penetration = amrex::min(1.0, amrex::max(d_xhi, d_yhi));

        const Real T = prims(i,j,k,cls.QT); // dT =0
        const Real rho = prims(i,j,k,cls.QRHO);
        Real Y[NUM_SPECIES] = {1.0};

        // Convert the target pressure to a target density at fixed T. This is
        // clearer and remains valid if the EOS is changed from ideal gas.
        Real target_rho = 0.0;
        cls.PYT2R(prob_parm.p_0, Y, T, target_rho);
        const Real drhodt = penetration * inverse_relaxation_time *(target_rho - rho);
        const Real u = prims(i,j,k,cls.QU);
        const Real v = prims(i,j,k,cls.QV);
        const Real w = prims(i,j,k,cls.QW);
        const Real total_specific_energy = prims(i,j,k,cls.QEINT) + 0.5*(u*u + v*v + w*w);

        rhs(i,j,k,cls.URHO) += drhodt;     
        rhs(i,j,k,cls.UMX) += u * drhodt;
        rhs(i,j,k,cls.UMY) += v * drhodt;
        rhs(i,j,k,cls.UMZ) += w * drhodt;
        rhs(i,j,k,cls.UET) += total_specific_energy * drhodt;
      }
    });

  };
};

} // namespace PROB

#endif
