#ifndef CNS_PROB_H_
#define CNS_PROB_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_ParmParse.H>

#include <PelePhysics.H>
#include <ReactorBase.H>

#include "Closures.h"
#include "RHS.h"

#if CNS_USE_EB    
#include <ebm.h>
#include <walltypes.h>
#endif

#include <bc_types.h>
#include <numbers>
#include <cmath>

namespace PROB {

// // LES closures
// struct LESparm {
//   // Smagorinsky constant
//   static constexpr Real Cs = 0.17; // Smag original 0.17 (off)
//   static constexpr int order = 2; // order of the numerical scheme for LES
//   static constexpr Real Scsgs = 0.4 ; //0.4 // turbulent Schmidt number
//   static constexpr Real Pr_o_Prsgs = 1.0; //0.1 // turbulent Prandtl number  
//   static constexpr bool fixDelta = false; // use fixed filter width
// };
// // Closure Index+Thermodynamics + Transport + LES + ATF
// typedef closures_dt<indicies_stat_t, transport_Pele_t , multispecies_pele_gas_t<indicies_t>,
//                     WALE_t<LESparm,indicies_t>, TFM_t<LESparm,indicies_t>>ProbClosures;

// ILES
using ProbClosures = closures_dt<indicies_t, transport_Pele_t , multispecies_pele_gas_t<indicies_t>>;

// problem parameters

struct ProbParm
{
  amrex::Real p0 = 590000; // total pressure [Pa]
  amrex::Real T0 = 1080;   // total temperature [K]
  amrex::Real M = 2.2;     // Mach number

  amrex::Real rho;
  amrex::Real u;
  amrex::Real ei;
  amrex::Real T;
  amrex::GpuArray<amrex::Real, NUM_SPECIES> Y = {0.0};

  amrex::Real rho_j;
  amrex::Real v_j;
  amrex::Real ei_j;
  amrex::Real T_j;
  amrex::GpuArray<amrex::Real, NUM_SPECIES> Y_jet = {0.0};

  amrex::Real cv_Tinf;

  bool spark = false; // add a hot spot in cavity floor (in bcnormal)

  bool record_statistics = false;    // record time avg (in prob_post_timestep)
  bool clean_aux_on_restart = false; // reset time avg (in prob_post_restart)

  bool do_bl = false;              // add BL at inflow (in bcnormal and prob_init)
  bool make_bl_on_restart = false; // set states before jet to be BL inflow condition (in prob_post_restart)

  bool make_init_on_restart = false; // reset states before x_reset to be inflow condition (in prob_post_restart)
  amrex::Real x_reset; 

  ProbParm () {
    // Inflow conditions
    Y[H2O_ID] = 1.068e-7 * T0 * T0 - 6.72e-5 * T0 + 2.986e-2; // fit from experiment
    Y[O2_ID] = 4.0 / 15.0 * (1.0 - Y[H2O_ID]);
    Y[N2_ID] = 11.0 / 15.0 * (1.0 - Y[H2O_ID]);

    // M > 1, terate to find gamma
    amrex::Real gamma = 1.313, p;
    auto eos = pele::physics::PhysicsType::eos();
    for (int iter = 0; iter < 10; ++iter) {
      // Isentropic relations
      T = T0 / (1 + 0.5 * (gamma - 1.0) * M * M);
      p = p0 * std::pow(1 + 0.5 * (gamma - 1) * M * M, -gamma / (gamma - 1));
      eos.PYT2R(p * pres_si2cgs, Y.begin(), T, rho);
      eos.RTY2G(rho, T, Y.begin(), gamma);
    }

    amrex::Real cv;
    eos.TY2Cv(T, Y.begin(), cv);
    cv_Tinf = cv * T * specenergy_cgs2si;
    eos.PYT2RE(p, Y.begin(), T, rho, ei);
    rho *= rho_cgs2si;
    ei *= specenergy_cgs2si;
    amrex::Real cs;
    eos.RTY2Cs(rho, T, Y.begin(), cs);
    u = M * cs * speed_cgs2si;

    amrex::Print() << "Inflow (gamma, rho, T, p, cv) = " << gamma << ", " << rho << ", "
                   << T << ", " << p << ", " << cv_Tinf / T << '\n';

    // Fuel conditions
    Y_jet[H2_ID] = 1.0;

    p0 = 845.0e4 + (755.e4 - 845.e4) * (T0 - 1100.0) / 300.0; // linearly varying
    M = 1.0;
    T0 = 288.0;
    gamma = 1.405;
    for (int iter = 0; iter < 10; ++iter) {
      T_j = T0 / (1 + 0.5 * (gamma - 1));
      p = p0 * std::pow(1 + 0.5 * (gamma - 1), -gamma / (gamma - 1));
      eos.PYT2R(p * pres_si2cgs, Y_jet.begin(), T_j, rho_j);
      eos.RTY2G(rho_j, T_j, Y_jet.begin(), gamma);
    }

    eos.PYT2RE(p, Y_jet.begin(), T_j, rho_j, ei_j);
    rho_j *= rho_cgs2si;
    ei_j *= specenergy_cgs2si;
    eos.RTY2Cs(rho_j, T_j, Y_jet.begin(), cs);
    v_j = M * cs * speed_cgs2si;

    amrex::Print() << "Fuel (gamma, rho, T, p) = " << gamma << ", " << rho_j << ", "
                   << T_j << ", " << p << '\n';
  }
};

// spark parametrs
// struct SparkParm{  
//   const Real t0 = 0.0085;          // spark time 
//   const Real x0 = 0.0;          // spark position
//   const Real y0 = 0.0;
//   const Real z0 = 0.08;     
//   const Real a  = 4.0*std::sqrt(std::log(10));
//   const Real Tmax    = 3000.0;
//   const Real energy  = 2000.0e-3; // 100 mJ
//   const Real pi = 3.14159265359;
//   //const Real ds   = std::sqrt(a/pi)*(energy);
//   const Real ds    = 4.0e-3;    // 5  mm
//   const Real dt    = 1.0e-3;    // 0.5 ms
//   const Real dt2   = dt*dt;     // 
//   const Real ds2   = ds*ds;     // 
//   const Real o_Volt = 0.25/(pi*pi*ds*ds*ds*dt);
//   const Real Cp0    =  1224; // [J/(kg K) assumed room temperature and phi=0.5
//   const Real  T0    = 300;
//   const Real  rho0  = 0.97 ;
//   const Real  ds0   = sqrt(a/pi)*std::pow(energy/(rho0*Cp0*(Tmax-T0)), 1.0/3.0);
//   // ds size that corresponds to a maximum temperature of Tmax, with given energy
// };


// numerical method parameters 
struct skewparm_t {
  static constexpr bool dissipation = true;         // no dissipation
  static constexpr int  order = 4;                  // order numerical scheme   (2 or 4)
  static constexpr amrex::Real C2skew=0.5,C4skew=0.016;    // Skew symmetric default  (0.1 0.016)
};

// struct viscous_param_t {
//   static constexpr int order = 2;                  // order numerical scheme   
//   static constexpr bool use_LES= false;
// };

struct wall_param {
  // static constexpr Real Twall = 285.5;                // wall temperature (if isothermal used)
  static constexpr bool solve_diffwall = true;      // solve viscous effects at walls
};

template <typename cls_t> class user_source_t;

// USED
typedef rhs_dt<skew_t<skewparm_t, ProbClosures>,
               viscousLES_t<user_source_t<ProbClosures>, ProbClosures>,
               reactor_sourceLES_t<user_source_t<ProbClosures>, ProbClosures>> ProbRHS;
// typedef rhs_dt<weno_t<ReconScheme::Teno5, ProbClosures>, viscous_t<user_source_t<ProbClosures>, ProbClosures>, reactor_source_t<user_source_t<ProbClosures>,ProbClosures >> ProbRHS;
// typedef rhs_dt<skew_t<skewparm_t, ProbClosures>, viscous_t<user_source_t<ProbClosures>, ProbClosures>, reactor_source_t<user_source_t<ProbClosures>,ProbClosures >> ProbRHS;

// define type of wall and EBM class
typedef adiabatic_wall_t<ProbClosures> TypeWall;
typedef ebm_t<TypeWall,wall_param,ProbClosures> ProbEB;
typedef manual_bc_t<ProbClosures> GlobalBC;

void inline inputs() {
  //	
}

// initial condition
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void prob_initdata(
    int i, int j, int k, Array4<Real> const& state,
    GeometryData const& geomdata, ProbClosures const& cls, ProbParm const& pp) {
  // Geometry
  const amrex::Real* prob_lo = geomdata.ProbLo();
  const amrex::Real* dx = geomdata.CellSize();
  // const amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];
  const amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1];
  // const amrex::Real z = prob_lo[2] + (k + 0.5) * dx[2];

  // Ventilated air
  const amrex::Real rho = pp.rho;
  const amrex::Real u = pp.u * (y > 0.0);
  const amrex::Real ei = pp.ei;
  state(i, j, k, cls.UMX) = rho * u;
  state(i, j, k, cls.UMY) = 0.0;
  state(i, j, k, cls.UMZ) = 0.0;
  state(i, j, k, cls.UET) = rho * ei + 0.5 * rho * u * u;
  for (int n = 0; n < NUM_SPECIES; ++n)
    state(i, j, k, cls.UFS + n) = rho * pp.Y[n];
}

/////////////////////////////// BC /////////////////////////////////////////////
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void bcnormal(
    const amrex::Real x[AMREX_SPACEDIM], amrex::Real dratio,
    const amrex::Real s_int[ProbClosures::NCONS],
    const amrex::Real s_refl[ProbClosures::NCONS],
    amrex::Real s_ext[ProbClosures::NCONS], const int idir, const int sgn,
    const amrex::Real time, GeometryData const& /*geomdata*/,
    ProbClosures const& cls, ProbParm const& pp) {
  amrex::Real rho, u, v, w, ei, Y[NUM_SPECIES];

  if (idir == 0) {
    // xlo: inflow
    rho = pp.rho;
    u = pp.u;
    v = 0.0;
    w = 0.0;
    ei = pp.ei;
    for (int n = 0; n < NUM_SPECIES; ++n) Y[n] = pp.Y[n];
  } else if (idir == 1) {
    // ylo
    if (x[0] * x[0] + x[2] * x[2] <= 0.001245 * 0.001245) {
      // fuel jet
      rho = pp.rho_j;
      u = 0.0;
      v = pp.v_j;
      w = 0.0;
      ei = pp.ei_j;
      for (int n = 0; n < NUM_SPECIES; ++n) Y[n] = pp.Y_jet[n];
    } else {
      // adiabatic wall
      rho = 0.0;
      for (int n = 0; n < NUM_SPECIES; ++n) rho += s_refl[cls.UFS + n];
      u = -s_refl[cls.UMX] / rho;
      v = -s_refl[cls.UMY] / rho;
      w = -s_refl[cls.UMZ] / rho;
      ei = s_refl[cls.UET] / rho - 0.5 * (u * u + v * v + w * w);
      for (int n = 0; n < NUM_SPECIES; ++n) Y[n] = s_refl[cls.UFS + n] / rho;
    }
  } else {
    return;
  }

  s_ext[cls.UMX] = rho * u;
  s_ext[cls.UMY] = rho * v;
  s_ext[cls.UMZ] = rho * w;
  s_ext[cls.UET] = rho * ei + 0.5 * rho * (u * u + v * v + w * w);
  for (int n = 0; n < NUM_SPECIES; ++n) s_ext[cls.UFS + n] = rho * Y[n];
}

///////////////////////////////AMR//////////////////////////////////////////////
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void user_tagging(
    int i, int j, int k, int nt_level, auto& tagarr, const auto& sdatafab,
    const auto& ebflag, const auto& geomdata, const ProbParm& prob_parm,
    int level) {
  // Geometry
  const amrex::Real* prob_lo = geomdata.ProbLo();
  const amrex::Real* dx = geomdata.CellSize();
  const amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];
  const amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1];
  const amrex::Real z = prob_lo[2] + (k + 0.5) * dx[2];
  
  // Tag walls
  amrex::Real d_wall = 0.0254_rt - y;  // dist from upper wall only
  d_wall = amrex::min(d_wall, 0.01905_rt - std::abs(z));  // dist from left and right walls
  if (d_wall > 0.0_rt && d_wall < 0.002_rt && x < 0.13_rt) {
    tagarr(i, j, k) = amrex::TagBox::SET;
  }

  // Tag EB cuts
  if (ebflag(i, j, k).isSingleValued()) {
    tagarr(i, j, k) = amrex::TagBox::SET;
  }

  // Tag density and YH2 gradient
  const int UFS = ProbClosures::UFS;
  amrex::Real rho[6] = {0.0};
  for (int n = 0; n < NUM_SPECIES; ++n) {
    rho[0] += sdatafab(i + 1, j, k, UFS + n);
    rho[1] += sdatafab(i - 1, j, k, UFS + n);
    rho[2] += sdatafab(i, j + 1, k, UFS + n);
    rho[3] += sdatafab(i, j - 1, k, UFS + n);
    rho[4] += sdatafab(i, j, k + 1, UFS + n);
    rho[5] += sdatafab(i, j, k - 1, UFS + n);
  }
  amrex::Real drhox = std::abs(rho[0] - rho[1]);
  amrex::Real drhoy = std::abs(rho[2] - rho[3]);
  amrex::Real drhoz = std::abs(rho[4] - rho[5]);
  if (amrex::max(drhox, drhoy, drhoz) > 0.1) {
    tagarr(i, j, k) = amrex::TagBox::SET;
  }

  amrex::Real YH2[6];
  YH2[0] = sdatafab(i + 1, j, k, UFS + H2_ID) / rho[0];
  YH2[1] = sdatafab(i - 1, j, k, UFS + H2_ID) / rho[1];
  YH2[2] = sdatafab(i, j + 1, k, UFS + H2_ID) / rho[2];
  YH2[3] = sdatafab(i, j - 1, k, UFS + H2_ID) / rho[3];
  YH2[4] = sdatafab(i, j, k + 1, UFS + H2_ID) / rho[4];
  YH2[5] = sdatafab(i, j, k - 1, UFS + H2_ID) / rho[5];
  amrex::Real dYH2x = std::abs(YH2[0] - YH2[1]);
  amrex::Real dYH2y = std::abs(YH2[2] - YH2[3]);
  amrex::Real dYH2z = std::abs(YH2[4] - YH2[5]);
  if (amrex::max(dYH2x, dYH2y, dYH2z) > 0.2) {
    tagarr(i, j, k) = amrex::TagBox::SET;
  }
}

///////////////////////////////SOURCE TERM /////////////////////////////////////
template <typename cls_t>
class user_source_t {
  public:

  // ATF options
  bool static constexpr use_ATF = false; // use adaptive thickening factor
  static constexpr int ATF_model = 1;   // 1: Classic  Colin/Charlette , 2: Rathore transformation
  //

  //...
  
  // viscous LES options
  static constexpr bool use_LES = false;         // use LES model in viscous term
  static constexpr int  order = 2;               // order numerical scheme   

  // compute chemistry
  bool static constexpr do_reactions = true;  // <<<<<<<<<<<<<<<<<<<<<<<<<
  bool static constexpr mask_cells_boundary = true; // avoid compute reactions in cells partially covered
  //

  

  // to use as a user source term, the function name must be src:
  // void inline src(const Geometry& geomdata, const amrex::MFIter &mfi,
  //                 const amrex::Array4<const amrex::Real> &prims,
  //                 const amrex::Array4<amrex::Real> &rhs, const cls_t *cls_d,
  //                 amrex::Real dt){

 
  // to use as a user source term after calling reactor, the function name must be rsrc:
  void inline rsrc(const Geometry& geomdata, const amrex::MFIter &mfi,
                  const amrex::Array4<const amrex::Real> &prims,
                  const amrex::Array4<amrex::Real> &rhs, const cls_t *cls_d,
                  amrex::Real dt, amrex::Real real_time){

  //   const Box& bxg = mfi.tilebox();
    
  //   // use device-friendly arrays instead of pointers
  //   //  const Real *prob_lo = geomdata.ProbLo();
  //   //  const Real *dx = geomdata.CellSize();
  //   auto prob_lo = geomdata.ProbLoArray();
  //   auto dx     = geomdata.CellSizeArray();
                
  //   // NTNUglobalnumbers const combustor;

  //   const Real tau_relax = 5.e-6;
  //   const Real coef =dt/tau_relax;
 
  //   const Real zbuffer  = 0.16;
  //   // const Real p_0      = combustor.pressure;

  //   // ------------------------------------------------------------------------------------
  //   amrex::ParallelFor(bxg,
  //     [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
  //     {
      
  //     // coordinates 
  //     // const Real x = prob_lo[0] + (i + Real(0.5)) * dx[0];
  //     // const Real y = prob_lo[1] + (j + Real(0.5)) * dx[1];
  //     const Real z = prob_lo[2] + (k + Real(0.5)) * dx[2];
  //     // const Real r = sqrt(x*x + y*y + (z-zexit)*(z-zexit));

  //     // get values from primitives
  //     const Real pres = prims(i,j,k,cls_t::QPRES);
  //     const Real T    = prims(i,j,k,cls_t::QT   ); 
  //     const Real rho  = prims(i,j,k,cls_t::QRHO );     
  //     Real Y[NUM_SPECIES]; for (int sp = 0; sp < NUM_SPECIES; sp ++) { Y[sp] = prims(i, j, k, cls_t::QFS + sp);}
  //     // u v w      
  //     const Real u  = prims(i,j,k,cls_t::QU );
  //     const Real v  = prims(i,j,k,cls_t::QV );
  //     const Real w  = prims(i,j,k,cls_t::QW );
  //     // Kinetic and Total energy (per density unit)  [J/rho] 
  //     const Real Et  = prims(i,j,k,cls_t::QEINT) + 0.5*(u*u + v*v + w*w);

  //     // pressure relaxation towards P0   
  //     const Real dP = (p_0- pres)*coef; Real drho;
  //     cls_d->PYT2R(dP,Y,T,drho);      // associated density change
  //     const Real drhodt = drho/dt;  // rate of change
   
  //     // buffer zone where relaxation applies 
  //     const bool buffer = (z > zbuffer); //&& (r > 0.06);

  //     if (buffer){        
 
  //       rhs(i,j,k,cls_t::UMX) += u*drhodt;
  //       rhs(i,j,k,cls_t::UMY) += v*drhodt; 
  //       rhs(i,j,k,cls_t::UMZ) += w*drhodt;
  //       rhs(i,j,k,cls_t::UET) += Et*drhodt;
  //       for (int sp = 0; sp < NUM_SPECIES; sp++)  
  //         rhs(i,j,k, cls_t::UFS + sp) += Y[sp] * drhodt;  
  //     }                        

  //   });

  }
};

} // namespace PROB

#endif
