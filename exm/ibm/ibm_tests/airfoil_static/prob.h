#ifndef CNS_PROB_H_
#define CNS_PROB_H_

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_ParmParse.H>
#include <Closures.h>
#include <RHS.h>
#include <ibm_solver.h>
#include <Constants.h>
#include <ibm_walltypes.h>

using namespace amrex;
using namespace universal_constants;

namespace PROB {

// constants
static constexpr Real Mach     = 2.0;
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
  Real p_oo    = 101325.0;
  Real T_oo    = 300.0;
  Real rho_oo  = p_oo/(Rgas*T_oo);
  Real c_oo    = sqrt(gam*Rgas*T_oo);
  Real u_oo    = c_oo*Mach;
  Real eint_oo = rho_oo*Cv*T_oo;
  Real kin_oo  = 0.5*rho_oo*u_oo*u_oo;
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
  static constexpr Real C2skew = 1.5, C4skew = 0.016;
};

struct ibmparm_t {
  public:
  static constexpr int  interp_order = 1;
  static constexpr int  extrap_order = 1;
  static constexpr Real alpha = 0.4;

  static constexpr int  interp_order_surf = 1;
  static constexpr int  extrap_order_surf = 1;
  static constexpr Real alpha_surf = 0.4;

  static constexpr int  ghost_layers = 1;
  static constexpr bool interior_is_solid = true;
};

// CLOSURES
typedef closures_dt<indicies_t, transport_const_t<methodparm_t>,
                    calorifically_perfect_gas_t<indicies_t>> ProbClosures;

typedef rhs_dt<weno_t<ReconScheme::WenoZ5, ProbClosures>,
               viscous_t<methodparm_t, ProbClosures>, no_source_t> ProbRHS;

typedef ibm_adiabatic_noslip_wall_t<ibmparm_t,ProbClosures> TypeWall;
typedef ibm_solver_t<TypeWall,ibmparm_t,ProbClosures> ProbIB;

// Static geometry: no update needed
inline void update_geometry(Real /*time*/,
                            Vector<GeomType>& /*geom_a*/,
                            int /*ngeom*/) {}

void inline inputs() {
  amrex::Print() << " ****** Static Diamond Wedge Airfoil (Ma=2) ******* \n";
}

//////////////////////////// Initial conditions ////////////////////////////////
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void prob_initdata (int i, int j, int k, amrex::Array4<amrex::Real> const& state,
      amrex::GeometryData const& geomdata, ProbClosures const& cls, ProbParm const& pp) {

  const Real* prob_lo = geomdata.ProbLo();
  const Real* dx      = geomdata.CellSize();
  Real x = prob_lo[0] + (i+0.5_rt)*dx[0];

  // Uniform freestream initialization
  state(i, j, k, cls.URHO) = pp.rho_oo;
  state(i, j, k, cls.UMX)  = pp.rho_oo * pp.u_oo;
  state(i, j, k, cls.UMY)  = 0.0;
#if (AMREX_SPACEDIM == 3)
  state(i, j, k, cls.UMZ)  = 0.0;
#endif
  state(i, j, k, cls.UET)  = pp.rho_oo * pp.eint_oo + pp.kin_oo;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void user_tagging(int i, int j, int k, int nt, auto& tagfab, const auto &sdatafab,
                  const Array4<const unsigned char>& ibfab, const auto& geomdata,
                  const ProbParm& pp, int level) {

  // Refine close to body (at all levels up to 2)
  if (level < 3) {
    bool foundGP = false;
    for (int ii = -1; ii <= 1; ii++) {
      for (int jj = -1; jj <= 1; jj++) {
#if (AMREX_SPACEDIM == 3)
        for (int kk = -1; kk <= 1; kk++) {
          foundGP = ibfab(i+ii,j+jj,k+kk,1) || foundGP;
        }
#else
        foundGP = ibfab(i+ii,j+jj,k,1) || foundGP;
#endif
      }
    }
    tagfab(i,j,k) = foundGP;
  }
}

//////////////////////////// Boundary conditions ///////////////////////////////
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
bcnormal(const amrex::Real x[AMREX_SPACEDIM], amrex::Real dratio, const amrex::Real s_int[ProbClosures::NCONS],
         const amrex::Real s_refl[ProbClosures::NCONS], amrex::Real s_ext[ProbClosures::NCONS],
         const int idir, const int sgn, const amrex::Real time,
         amrex::GeometryData const& /*geomdata*/,  ProbClosures const& closures, ProbParm const& pp)
{
  const int URHO = ProbClosures::URHO;
  const int UMX  = ProbClosures::UMX;
  const int UMY  = ProbClosures::UMY;
  const int UET  = ProbClosures::UET;
  const int face = (idir+1)*sgn;

  switch(face)
  {
    case  1:  // x-low: supersonic inflow
      s_ext[URHO] = pp.rho_oo;
      s_ext[UMX]  = pp.rho_oo * pp.u_oo;
      s_ext[UMY]  = 0.0;
#if (AMREX_SPACEDIM == 3)
      s_ext[ProbClosures::UMZ] = 0.0;
#endif
      s_ext[UET]  = pp.rho_oo * pp.eint_oo + pp.kin_oo;
      break;
    default:  // outflow
      break;
  }
}

}
#endif
