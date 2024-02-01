#ifndef CNS_PROB_H_
#define CNS_PROB_H_

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_ParmParse.H>
#include <Closures.h>

using namespace amrex;

namespace PROB {

////////////////////////////////EQUATIONS///////////////////////////////////////
constexpr int do_pde = 0; // 0=CNS, 1=CNS-NLDE
// Select the variables to solve/store and write
// Independent (solved) variables
#define URHO  0
#define UMX   1
#define UMY   2
#define UMZ   3
#define UET   4
#define NCONS 5

// Dependent (derived) variables
#define QRHO   0
#define QU     1
#define QV     2
#define QW     3
#define QT     4
#define QPRES  5
#define NPRIM  6

// fluxes
// typedef euler_derived_base_t<visc_suth_t, cond_suth_t, calorifically_perfect_gas_t> ProbClosures;

// typedef 
////////////////////////////////////////////////////////////////////////////////

//////////////////////////////DISCRETISATION////////////////////////////////////
#define NGHOST 3 // TODO: make this an automatic parameter?
void inline inputs() {

  ParmParse pp;

  // Numerical operators
  //-1 = N/A (Incase of periodic)
  // 0 = Interior           3 = Symmetry
  // 1 = Inflow             4 = SlipWall
  // 2 = Outflow            5 = NoSlipWall
  // 6 = user defined
  pp.addarr("cns.lo_bc", std::vector<int>{2,-1,-1});
  pp.addarr("cns.hi_bc", std::vector<int>{2,-1,-1});
  pp.add   ("cns.order_rk", 3); // -2, 1, 2 or 3"
  pp.add   ("cns.stages_rk", 3); // 1, 2 or 3
  pp.add   ("cns.rhs_euler", 1); // 0=false, 1=true
  pp.add   ("cns.rhs_visc", 0); // 0=false, 1=true
  pp.add   ("cns.rhs_source", 0); // 0=false, 1=true
  pp.add   ("cns.flux_euler", 2); // 0=riemann solver, 1=KEEP/AD, 2=WENO5
  pp.add   ("cns.order_keep", 4); // Order of accuracy=2, 4 or 6"
  pp.add   ("cns.art_diss", 0); // 0=none, 1=artificial dissipation
  // pp.add   ("cns.nghost",2);
  pp.add   ("cns.screen_output", 1); // 0=quiet, 1=verbose
  pp.add   ("cns.verbose", 1); // 0=quiet, 1=verbose

  // lets add viscosity type and different flux dissipations!


  // debugging
  // pp.add("amrex.fpe_trap_invalid",1);
  // pp.add("amrex.fpe_trap_zero",1);
  // pp.add("amrex.fpe_trap_overflow",1);

}
////////////////////////////////////////////////////////////////////////////////

///////////////////////////////CLOSURES/////////////////////////////////////////
typedef closures_derived_base_t<visc_suth_t, cond_suth_t, calorifically_perfect_gas_t> ProbClosures;
// user can also define their own closure class and use it here by naming it ProbClosures
// template <typename Visc, typename Cond, typename Thermo>
// class closures_derived_user_t : public Cond, public Visc, public Thermo
// {
  // private:
  //
  // public:
// };

// problem parameters
struct ProbParm
{
  Real p_l = 1.0;
  Real p_r = 0.1;
  Real rho_l = 1.0;
  Real rho_r = 0.125;
  Real u_l = 0.0;
  Real u_r = 0.0;
};

// initial condition
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void prob_initdata (int i, int j, int k, Array4<Real> const& state, GeometryData const& geomdata, ProbClosures const& closures, ProbParm const& prob_parm) {
  const Real* prob_lo = geomdata.ProbLo();
  const Real* prob_hi = geomdata.ProbHi();
  const Real* dx      = geomdata.CellSize();

  Real x = prob_lo[0] + (i+Real(0.5))*dx[0];
  Real Pt, rhot, uxt;
  if (x < prob_hi[0]/2) {
      Pt = prob_parm.p_l;
      rhot = prob_parm.rho_l;
      uxt = prob_parm.u_l;
  } else {
      Pt = prob_parm.p_r;
      rhot = prob_parm.rho_r;
      uxt = prob_parm.u_r;
  }
  state(i,j,k,URHO ) = rhot;
  state(i,j,k,UMX  ) = rhot*uxt;
  state(i,j,k,UMY  ) = Real(0.0);
  state(i,j,k,UMZ  ) = Real(0.0);
  Real et = Pt/(closures.gamma-Real(1.0));
  state(i,j,k,UET) = et + Real(0.5)*rhot*uxt*uxt;
}

// boundary conditions
/**
 * \brief Fill external boundary conditions for ghost cells.
 *
 * @param x         ghost cell cooridinates.
 * @param dr        wall-ghost/wall-first internal distance ratio 
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
AMREX_GPU_DEVICE AMREX_FORCE_INLINE 
void bcnormal(const Real x[AMREX_SPACEDIM], Real dratio, const Real s_int[NCONS],
         const Real s_refl[NCONS], Real s_ext[NCONS],
         const int idir, const int sgn, const Real time,
         GeometryData const& /*geomdata*/,  ProbClosures const& closures, ProbParm const& prob_parm)
{
  if (idir == 1) { // ylo or yhi

    Abort("bcnormal not coded");

    // Real q_ext[NPRIM] = {0.0};
    // no-slip
    // q_ext[QU]    = -s_int[UMX]/s_int[URHO];
    // q_ext[QV]    = -s_int[UMY]/s_int[URHO];
    // q_ext[QW]    = -s_int[UMZ]/s_int[URHO];

    // // dp/dn = 0
    // Real eint_int = (s_int[UET] - 0.5*(s_int[UMX]*s_int[UMX] + s_int[UMY]*s_int[UMY] + s_int[UMZ]*s_int[UMZ])/s_int[URHO])/s_int[URHO];
    // Real p_int = (parm.eos_gamma - 1.0)*s_int[URHO]*eint_int;
    // q_ext[QPRES] = p_int;
    // // T=Twall
    // Real T_int = p_int/(parm.Rspec*s_int[URHO]); 
    // q_ext[QT]    = max(prob_parm.Tw  +  dratio*(prob_parm.Tw - T_int),50.0);
    // // rho = eos(P,T)
    // q_ext[QRHO]  = q_ext[QPRES]/(parm.Rspec*q_ext[QT]);

    // // convert prims to cons
    // s_ext[URHO] = q_ext[QRHO];
    // s_ext[UMX] = q_ext[QRHO]*q_ext[QU];
    // s_ext[UMY] = q_ext[QRHO]*q_ext[QV];
    // s_ext[UMZ] = q_ext[QRHO]*q_ext[QW];
    // Real ekin_ext = 0.5*(q_ext[QU]*q_ext[QU] + q_ext[QV]*q_ext[QV] + q_ext[QW]*q_ext[QW]); 
    // Real eint_ext = q_ext[QPRES]/(q_ext[QRHO]*(parm.eos_gamma - 1.0));
    // s_ext[UET] = q_ext[QRHO]*(eint_ext + ekin_ext);
  }
}

// source term
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void user_source(int i, int j, int k, const auto& state, const auto& rhs, const ProbParm& lprobparm, ProbClosures const& closures, auto const dx) {

}
////////////////////////////////////////////////////////////////////////////////


///////////////////////////////AMR//////////////////////////////////////////////
AMREX_GPU_DEVICE AMREX_FORCE_INLINE 
void user_tagging(int i, int j, int k, int nt_level, auto& tagfab, const auto &sdatafab, const auto& geomdata, const ProbParm& prob_parm , int level) {

      Real dengrad_threshold = 0.5;
      Real drhox = Math::abs(sdatafab(i+1,j,k,URHO) - sdatafab(i-1,j,k,URHO))/sdatafab(i,j,k,URHO);
      if (drhox > dengrad_threshold) {
        tagfab(i,j,k) = true;
        tagfab(i+1,j,k) = true;
        tagfab(i+2,j,k) = true;
        tagfab(i+3,j,k) = true;
      }
  }
////////////////////////////////////////////////////////////////////////////////

} // namespace prob
#endif