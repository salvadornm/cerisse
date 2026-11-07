#ifndef CNS_PROB_H_
#define CNS_PROB_H_

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_ParmParse.H>
#include <AMReX_AmrLevel.H>
#include <Closures.h>
#include <RHS.h>
#include <ibm_solver.h>
#include <ratio>
#include <Constants.h>
#include <ibm_walltypes.h>

using namespace amrex;
using namespace universal_constants;

namespace PROB {

// constants
static constexpr Real Mach     = 4.0;
static constexpr Real Mw       = 28.96e-3;
static constexpr Real gam      = 1.4;
static constexpr Real Rgas     = gas_constant/Mw;
static constexpr Real Reynolds = 10000;
static constexpr Real Pr       = 0.7;
static constexpr Real Cv       = Rgas/(gam - 1.0);
static constexpr Real Cp       = gam*Cv;
static constexpr Real viscos   = 1.0/Reynolds;
static constexpr Real lambda   = viscos*Cp/Pr;

//////////////////////////// Physical modelling ////////////////////////////////
struct ProbParm
{
  // free-stream conditions
  Real p_oo    = 100000.0;
  Real T_oo    = 300.0;
  Real rho_oo  = p_oo/(Rgas*T_oo);
  Real c_oo    = sqrt(gam*Rgas*T_oo);
  Real u_oo    = c_oo*Mach;
  Real eint_oo = rho_oo*Cv*T_oo;
  Real kin_oo  = 0.5*rho_oo*u_oo*u_oo;

  // right state
  Real p_r     = p_oo;
  Real T_r     = T_oo;
  Real rho_r   = p_r/(Rgas*T_r);
  Real u_r     = 0.0;
  Real eint_r  = rho_r*Cv*T_r;

  // centre of geometry
  Real x0 = 0.0; Real y0 = 0.0;
#if (AMREX_SPACEDIM == 3)
  Real z0 = 0.0;
#endif
  // initial shock position
  //Real xshock = -0.9;
  Real xshock = -0.4;

};

struct methodparm_t {
  public:
  static constexpr int  order = 2;
  static constexpr Real conductivity = lambda;
  static constexpr Real viscosity    = viscos;
  static constexpr bool use_LES = false;
};

struct skewparm_t {
  public:
  static constexpr bool dissipation = true;
  static constexpr int  order = 4;
  static constexpr Real C2skew=1.5,C4skew=0.016;
};

struct ibmparm_t {
  public:
  static constexpr int  interp_order = 1;
  static constexpr int  extrap_order = 1;
  static constexpr Real alpha = 0.6;

  static constexpr int  interp_order_surf = 1;
  static constexpr int  extrap_order_surf = 1;
  static constexpr Real alpha_surf = 0.6;

  static constexpr int  ghost_layers = 1;
  static constexpr bool interior_is_solid = true;
};

// CLOSURES
typedef closures_dt<indicies_t, transport_const_t<methodparm_t>,
                    calorifically_perfect_gas_t<indicies_t>> ProbClosures;

typedef rhs_dt<weno_t<ReconScheme::WenoZ5, ProbClosures>, no_diffusive_t, no_source_t > ProbRHS;

typedef ibm_adiabatic_noslip_wall_t<ibmparm_t,ProbClosures> TypeWall;
typedef ibm_solver_t<TypeWall,ibmparm_t,ProbClosures> ProbIB;

// Static geometry: no update needed
inline void update_geometry(Real /*time*/,
                            Vector<GeomType>& /*geom_a*/,
                            int /*ngeom*/) {}

void inline inputs() {
  amrex::Print() << " ****** Phase 1 IBM Backend Verification Test ******* " << std::endl;
}

//////////////////////////// Initial conditions ////////////////////////////////
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void prob_initdata (int i, int j, int k, amrex::Array4<amrex::Real> const& state,
      amrex::GeometryData const& geomdata, ProbClosures const& cls, ProbParm const& pparm) {

  const Real* prob_lo = geomdata.ProbLo();
  const Real* dx      = geomdata.CellSize();

  Real x = prob_lo[0] + (i+0.5_rt)*dx[0];
  Real rhot,eint,u[3]={0.0};

  if (x < pparm.xshock) {
    rhot =  pparm.rho_oo;
    u[0] =  pparm.u_oo;
    eint =  pparm.eint_oo;
  }
  else {
    rhot =  pparm.rho_r;
    u[0] =  pparm.u_r;
    eint =  pparm.eint_r;
  }

  state(i, j, k, cls.URHO) = rhot;
  state(i, j, k, cls.UMX)  = rhot * u[0];
  state(i, j, k, cls.UMY)  = Real(0.0);
  state(i, j, k, cls.UMZ)  = Real(0.0);
  state(i, j, k, cls.UET)  = eint + Real(0.5) * rhot * u[0] * u[0] ;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void user_tagging(int i, int j, int k, int nt, auto& tagfab, const auto &sdatafab,
                  const Array4<const unsigned char>& ibfab, const auto& geomdata,
                  const ProbParm& pparm , int level) {

  // refine close to body (at all levels)
  if (ibfab(i,j,k,1)) {
    for (int ii = -1; ii <= 1; ii++) {
      for (int jj = -1; jj <= 1; jj++) {
#if (AMREX_SPACEDIM == 3)
        for (int kk = -1; kk <= 1; kk++) {
          tagfab(i+ii,j+jj,k+kk) = true;
        }
#else
        tagfab(i+ii,j+jj,k) = true;
#endif
      }
    }
  }
}

//////////////////////////// Boundary conditions ///////////////////////////////
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
bcnormal(const amrex::Real x[AMREX_SPACEDIM], amrex::Real dratio, const amrex::Real s_int[ProbClosures::NCONS],
         const amrex::Real s_refl[ProbClosures::NCONS], amrex::Real s_ext[ProbClosures::NCONS],
         const int idir, const int sgn, const amrex::Real time,
         amrex::GeometryData const& /*geomdata*/,  ProbClosures const& closures, ProbParm const& pparm)
{
  const int URHO = ProbClosures::URHO;
  const int UMX  = ProbClosures::UMX;
  const int UMY  = ProbClosures::UMY;
  const int UET  = ProbClosures::UET;
  const int face = (idir+1)*sgn;

  switch(face)
  {
    case   1:  // WEST - inflow
      s_ext[URHO] = pparm.rho_oo;
      s_ext[UMX]  = pparm.rho_oo * pparm.u_oo;
      s_ext[UMY]  = 0.0;
#if (AMREX_SPACEDIM == 3)
      s_ext[ProbClosures::UMZ]  = 0.0;
#endif
      s_ext[UET]  = pparm.eint_oo + pparm.kin_oo;
      break;
    default:
      break;
  }
}

}
#endif
