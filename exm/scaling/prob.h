#ifndef PROB_H
#define PROB_H

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_ParmParse.H>
#include <Closures.h>
#include <RHS.h>
#include <Constants.h>

using namespace amrex;

namespace PROB {

static constexpr Real Reynolds = 1600.0;
static constexpr Real Mach     = 0.1;    
static constexpr Real Prandtl  = 0.71;  
static constexpr Real gam      = 1.4;  
static constexpr Real Rgas     = gas_constant*1000.0/28.96;  
static constexpr Real Cp       = Rgas*gam/(gam - 1.0) ;  

struct viscparm_t {

  public :

  static constexpr int  order = 2;
  static constexpr bool use_LES = false;

  // constant viscosity and conductivity  (infinite Reynolds)
  static constexpr Real viscosity    = 1.0/Reynolds;  
  static constexpr Real conductivity = viscosity*Cp/Prandtl;

};

using ProbClosures = closures_dt<indicies_t, transport_const_t<viscparm_t>,
                          calorifically_perfect_gas_t<indicies_t>>;
//////////////////////////////DISCRETISATION////////////////////////////////////

// WENO/TENO (TENO5 works  Ma=0.1, TENO6 not working yet in this case)
using ProbRHS      =   rhs_dt<weno_t<ReconScheme::Teno5, ProbClosures>,
                            viscous_t<viscparm_t, ProbClosures>,no_source_t>;

// KEEP (works test case Ma=0.1)
// using ProbRHS      =   rhs_dt<keep_euler_t<false,false,4, ProbClosures>,
//                             viscous_t<viscparm_t, ProbClosures>,no_source_t>;

// SKEW
// struct skewparm_t {

//   public:

//   static constexpr bool dissipation = true;         // no dissipation
//   static constexpr int  order = 4;                  // order numerical scheme
//   static constexpr Real C2skew=0.05,C4skew=0.016;   // Skew symmetric default

// };

// using ProbRHS      =   rhs_dt<skew_t<skewparm_t, ProbClosures>,
//                             viscous_t<viscparm_t, ProbClosures>,no_source_t>;


// CD
//using ProbRHS      =   rhs_dt<centraldif_t<false,false,4, ProbClosures>,
//                            viscous_t<viscparm_t, ProbClosures>,no_source_t>;


// problem parameters
struct ProbParm {
  //bool convecting = false;

  Real omega_x = Real(1.0);   // [rad s^-1]
  Real omega_y = Real(1.0);   // [rad s^-1]
  Real omega_z = Real(1.0);   // [rad s^-1]
  Real L = Real(1.0);         // [m]
  Real T0 = Real(300.0);        // [K]
  Real c0 = std::sqrt(gam * Rgas*T0); // soud speed [m/s]
  Real v0 = Mach*c0;            // ref velocity
  Real mu0 = 1.0/Reynolds;
  Real rho0 = Reynolds*mu0/(v0*L);  // density [kg/m3] (from Reynolds) 
  Real time0 = L/v0;            // ref time [s]
  Real p0 = rho0*Rgas*T0;       // pressure [Pa]   
};
///////////////////////////////////////////////////////////////////////
void inline inputs() {
  ProbParm data;


  Real pref = data.rho0*data.v0*data.v0;

  amrex::Print() << "**************  " << std::endl;
  amrex::Print() << " Taylor Green  Test  " << std::endl;
  amrex::Print() << " Ma  =  " << Mach << std::endl;
  amrex::Print() << " Re (imposed)  =  " << Reynolds << std::endl;  
  amrex::Print() << " Re (calc)     =  " << data.rho0*data.v0*data.L/data.mu0 << std::endl;  
  amrex::Print() << " rho0 =  " << data.rho0 << " p0= " << data.p0 <<std::endl; 
  amrex::Print() << " v0 =    " << data.v0   << " T0= " << data.T0 <<std::endl; 
  amrex::Print() << " c0 (calc) =    " << sqrt(gam*data.p0/data.rho0)  << " c0= " << data.c0 <<std::endl; 
  amrex::Print() << " rho0 v0^2/p0 =  " << pref/data.p0 << std::endl;
  
  amrex::Print() << " time0 =    " << data.time0   << std::endl; 
  
  amrex::Print() << "**************  " << std::endl;
}

// initial condition
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
prob_initdata(int i, int j, int k, Array4<Real> const &state,
              GeometryData const &geomdata, ProbClosures const &cls,
              const ProbParm &pparm) {
  // Geometry
  const Real *prob_lo = geomdata.ProbLo();
  const Real *dx = geomdata.CellSize();
  const Real x = prob_lo[0] + (i + 0.5) * dx[0];
  const Real y = prob_lo[1] + (j + 0.5) * dx[1];
  const Real z = prob_lo[2] + (k + 0.5) * dx[2];

  // ref values
  Real v0   = pparm.v0;  
  Real rho0 = pparm.rho0;
  Real p0   = pparm.p0;

  // TGV functions
  Real u[3] = {0.0};
  u[0] = v0 * sin(pparm.omega_x * x / pparm.L) *
         cos(pparm.omega_y * y / pparm.L) * cos(pparm.omega_z * z / pparm.L);
  u[1] = -v0 * cos(pparm.omega_x * x / pparm.L) *
         sin(pparm.omega_y * y / pparm.L) * cos(pparm.omega_z * z / pparm.L);
  // if (pparm.convecting) {
  //   u[0] += v0;
  //   u[1] += v0;
  // }
  const Real p = p0 + rho0 * v0 * v0 / Real(16.0) *
                          (cos(2.0 * pparm.omega_x * x / pparm.L) +
                           cos(2.0 * pparm.omega_y * y / pparm.L)) *
                          (cos(2.0 * pparm.omega_z * z / pparm.L) + Real(2.0));
  Real rho = p / (cls.Rspec * pparm.T0);
  Real eint = cls.cv * pparm.T0;



  // Set the state
  state(i, j, k, ProbClosures::URHO) = rho;
  state(i, j, k, ProbClosures::UMX) = rho * u[0];
  state(i, j, k, ProbClosures::UMY) = rho * u[1];
  state(i, j, k, ProbClosures::UMZ) = rho * u[2];
  state(i, j, k, ProbClosures::UET) =
      rho * (eint + Real(0.5) * (u[0] * u[0] + u[1] * u[1] + u[2] * u[2]));
}

// boundary conditions
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
bcnormal(const Real x[AMREX_SPACEDIM], Real dratio, const Real s_int[ProbClosures::NCONS],
         const Real s_refl[ProbClosures::NCONS], Real s_ext[ProbClosures::NCONS], const int idir,
         const int sgn, const Real time, GeometryData const & /*geomdata*/,
         ProbClosures const &cls, ProbParm const &pparm) {
  Abort("bcnormal not set");
}

// source term
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
user_source(int i, int j, int k, const auto &state, const auto &rhs,
            const ProbParm &pparm, ProbClosures const &cls, auto const dx) {}
////////////////////////////////////////////////////////////////////////////////

///////////////////////////////AMR//////////////////////////////////////////////
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
user_tagging(int i, int j, int k, int nt_level, auto &tagfab,
             const auto &sdatafab, const auto &geomdata, const ProbParm &pparm,
             int level) {}
////////////////////////////////////////////////////////////////////////////////

} // namespace PROB
#endif
