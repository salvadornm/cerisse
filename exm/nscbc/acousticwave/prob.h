#ifndef CNS_PROB_H_
#define CNS_PROB_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_ParmParse.H>
#include <Closures.h>
#include <RHS.h>

using namespace amrex;

namespace PROB {

// problem parameters
struct ProbParm {

  Real rho0= 1.18;  // kg/m3    
  Real gamma = 1.4;
  Real L0 = 1;
  Real u0 = 1.0;
  Real p0 = 101325.0; // 1 atm
  Real c0 = std::sqrt(gamma * p0 / rho0);
  Real t0 = L0/(u0 + c0);  // characteristic time scale
  // wave parameters
  Real A = 5.0;
  Real B = 10.0;
  
  Real sigma = 0.1;

};

// numerical method parameters
struct methodparm_t {

  public:

  static constexpr bool dissipation = true;         // no dissipation
  static constexpr int  order = 6;                  // order numerical scheme
  static constexpr Real C2skew=0.1,C4skew=0.016;   // Skew symmetric default

};

typedef closures_dt<indicies_t, visc_suth_t, cond_suth_t,
                    calorifically_perfect_gas_t<indicies_t>>
    ProbClosures;

// HLLC Riemann solver    
//typedef rhs_dt<riemann_t<false, ProbClosures>, no_diffusive_t, no_source_t    > ProbRHS;
// skew-symmetric
typedef rhs_dt<skew_t<methodparm_t, ProbClosures>, no_diffusive_t, no_source_t > ProbRHS;
// KEEP 2/4/6
//typedef rhs_dt<keep_euler_t<false,false,4, ProbClosures>, no_diffusive_t,  no_source_t> ProbRHS;
// WENO & TENO   WenoZ5/Teno5/Teno6
//typedef rhs_dt<weno_t<ReconScheme::Teno5, ProbClosures>, no_diffusive_t,  no_source_t>  ProbRHS;
    
//
void inline inputs() {
  ProbParm data;
  amrex::Print() << "**************                          ************ " << std::endl;
  amrex::Print() << " One-dimensional acoustic wave propagation " << std::endl;
  amrex::Print() << "**************                          ************ " << std::endl;
  amrex::Print() << " t0 [s]  =  " << data.t0 << std::endl;
  Real Ma = data.u0 / data.c0;     
  amrex::Print() << " Mach   =  " << Ma << std::endl;
  amrex::Print() << "**************                          ************ " << std::endl;
}

// initial condition
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
prob_initdata(int i, int j, int k, Array4<Real> const &state,
              GeometryData const &geomdata, ProbClosures const &cls,
              ProbParm const &prob_parm) {
  const Real *prob_lo = geomdata.ProbLo();
  const Real *prob_hi = geomdata.ProbHi();
  const Real *dx = geomdata.CellSize();

  Real x = prob_lo[0] + (i + Real(0.5)) * dx[0];
  Real Pt, rhot, uxt,arg;
    

  Real c0 = prob_parm.c0;

  arg   = prob_parm.B*(x - Real(0.5)*prob_parm.L0)/prob_parm.L0;

  // acoustic wave
  Real dp = 1.0e-3 * prob_parm.p0 * exp(-arg*arg); // small
  Pt   = prob_parm.p0 + dp;
  rhot = prob_parm.rho0 + dp/(c0*c0);
  uxt  = prob_parm.u0 + dp/(prob_parm.rho0*c0);
  //

  state(i, j, k, cls.URHO) = rhot;
  state(i, j, k, cls.UMX) = rhot * uxt;
  state(i, j, k, cls.UMY) = Real(0.0);
  state(i, j, k, cls.UMZ) = Real(0.0);
  Real et = Pt / (cls.gamma - Real(1.0));
  state(i, j, k, cls.UET) = et + Real(0.5) * rhot * uxt * uxt;
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
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
bcnormal(const Real x[AMREX_SPACEDIM], Real dratio, const Real s_int[5],
         const Real s_refl[ProbClosures::NCONS], Real s_ext[5], const int idir,
         const int sgn, const Real time, GeometryData const & /*geomdata*/,
         ProbClosures const &closures, ProbParm const &prob_parm) {



}

// source term
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
user_source(int i, int j, int k, const auto &state, const auto &rhs,
            const ProbParm &lprobparm, ProbClosures const &closures,
            auto const dx) {}
////////////////////////////////////////////////////////////////////////////////

///////////////////////////////AMR//////////////////////////////////////////////
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
user_tagging(int i, int j, int k, int nt_level, auto &tagfab,
             const auto &sdatafab, const auto &geomdata,
             const ProbParm &prob_parm, int level) {

}
////////////////////////////////////////////////////////////////////////////////

} // namespace PROB
#endif
