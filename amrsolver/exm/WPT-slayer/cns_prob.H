#ifndef CNS_PROB_H_
#define CNS_PROB_H_

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_PROB_AMR_F.H>
#include <AMReX_ParmParse.H>
#include <Closures.h>
#include <AMReX_Print.H>

using namespace amrex;


#define URHO  0
#define UMX   1
#define UMY   2
#define UMZ   3
#define UET   4
#define NCONS 5

#define QRHO   0
#define QU     1
#define QV     2
#define QW     3
#define QT     4
#define QPRES  5
#define NPRIM  6

#define NGHOST 3
namespace PROB {

struct ProbParm
{
  Real M1 = 4;              // upper Mach number
  Real Re1 = 400;           // upper Reynolds number
  Real U1 = 1200;           // upper velocity (m/s)
  Real Uref = U1/M1;        // reference velocity (m/s) - upper speed of sound              
  Real U2 = 0.375*U1;       // lower velocity (m/s)
  Real T1 = 100;            // upper temperature (K)
  Real T2 = 3*T1;           // lower temperature (K)
  Real dudymax = 1e3;           // max mean velocity gradient (1/s)
  Real dw = (U1-U2)/dudymax;     // vorticity thickness (m)
  Real lref = dw/2;              // reference length (m)
  Real mu1 = (1.458e-6)*(pow(T1,(1.5)))/(110.4+T1);  // Sutherland's model of viscosity (Pa s)
  Real rho1 = Re1*mu1/(Uref*lref);
  Real P1 = rho1*287*T1;
    // sod problem parameters below
  Real p_l = 1.0;
  Real p_r = 0.1;
  Real rho_l = 1.0;
  Real rho_r = 0.125;
  Real u_l = 0.0;
  Real u_r = 0.0;
};

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
////////////////////////////////////////////////////////////////////////////////


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
  pp.add   ("cns.rhs_visc", 1); // 0=false, 1=true
  pp.add   ("cns.rhs_source", 0); // 0=false, 1=true
  pp.add   ("cns.flux_euler", 3); // 0=riemann solver, 1=KEEP/AD, 2=WENO5, 3=NLDE
  pp.add   ("cns.screen_output", 1); // 0=quiet, 1=verbose
  pp.add   ("cns.verbose", 1); // 0=quiet, 1=verbose

}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void prob_initdata (int i, int j, int k, amrex::Array4<amrex::Real> const& state, amrex::GeometryData const& geomdata, ProbClosures const& closures, ProbParm const& prob_parm) {
  const Real* prob_lo = geomdata.ProbLo();
  const Real* prob_hi = geomdata.ProbHi();
  const Real* dx      = geomdata.CellSize();

  // Print() << "U1=" << prob_parm.U1 << "\n"; // print value of U1 to screen ... why does this print more than just once?!
  // Print() << "U2=" << prob_parm.U2 << "\n";
  // Print() << "dudymax=" << prob_parm.dudymax << "\n";
  // Print() << "dw=" << prob_parm.dw << "\n";
  // Print() << "lref=" << prob_parm.lref << "\n";
  // Print() << "mu1=" << prob_parm.mu1 << "\n";
  // Print() << "T1=" << prob_parm.T1 << "\n";
  // Print() << "T2=" << prob_parm.T2 << "\n";
  // Print() << "rho1=" << prob_parm.rho1 << "\n";

  Real x = prob_lo[0] + (i+Real(0.5))*dx[0];
  Real y = prob_lo[1] + (j+Real(0.5))*dx[1];
  Real Pt, rhot, uxt, uyt, Tt, u1t, u2t, T1t, T2t;
  
  // mean pressure is uniform
  Pt = prob_parm.P1;

  // analytical expression for velocity distribution
  u1t = prob_parm.U1;
  u2t = prob_parm.U2;
  uxt = Real(0.5)*((u1t+u2t)+(u1t-u2t)*tanh(Real(2.0)*y/prob_parm.dw));
  
  // ananlytical expression for temperature distribution
  T1t = prob_parm.T1;
  T2t = prob_parm.T2;
  Tt = T1t*(T2t/T1t)*((1-(uxt/u1t))/(1-(u2t/u1t)))
      +T1t*(((uxt/u1t)-(u2t/u1t))/(1-(u2t/u1t)))
      +T1t*Real(0.5)*(closures.gamma-Real(1.0))*prob_parm.M1*prob_parm.M1*(1-(uxt/u1t))*((uxt/u1t)-(u2t/u1t));
  
  // density profile assuming perfect gas
  rhot = Pt/(closures.Rspec*Tt);

  state(i,j,k,URHO ) = rhot;
  state(i,j,k,UMX  ) = rhot*uxt;
  state(i,j,k,UMY  ) = Real(0.0);
  state(i,j,k,UMZ  ) = Real(0.0);
  Real et = Pt/(closures.gamma-Real(1.0));
  state(i,j,k,UET) = et + Real(0.5)*rhot*uxt*uxt;
  
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE 
void user_tagging(int i, int j, int k, int nt_lev, auto& tagfab, const auto &sdatafab, const auto& geomdata, const ProbParm& prob_parm , int level) {

// this section is for AMR refinement - fine as is.
      Real dengrad_threshold = 0.5;
      amrex::Real drhox = amrex::Math::abs(sdatafab(i+1,j,k,URHO) - sdatafab(i-1,j,k,URHO))/sdatafab(i,j,k,URHO);
      if (drhox > dengrad_threshold) {
        tagfab(i,j,k) = true;
        tagfab(i+1,j,k) = true;
        tagfab(i+2,j,k) = true;
        tagfab(i+3,j,k) = true;
      }
  }


AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void user_source(int i, int j, int k, const auto& state, const auto& rhs, const ProbParm& lprobparm, ProbClosures const& closures, auto const dx) {

}

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
// this section for setting wall bc ghost cells. not needed for shear layer prob, but will later for blayer.
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
bcnormal(const amrex::Real x[AMREX_SPACEDIM], amrex::Real dratio, const amrex::Real s_int[NCONS],
         const amrex::Real s_refl[NCONS], amrex::Real s_ext[NCONS],
         const int idir, const int sgn, const amrex::Real time,
         amrex::GeometryData const& /*geomdata*/,  ProbClosures const& closures, ProbParm const& prob_parm)
{
  if (idir == 1) { // ylo or yhi

    amrex::Abort("bcnormal not coded");

    // Real q_ext[NPRIM] = {0.0};
    // no-slip
    // q_ext[QU]    = -s_int[UMX]/s_int[URHO];
    // q_ext[QV]    = -s_int[UMY]/s_int[URHO];
    // q_ext[QW]    = -s_int[UMZ]/s_int[URHO];

    // // dp/dn = 0
    // amrex::Real eint_int = (s_int[UET] - 0.5*(s_int[UMX]*s_int[UMX] + s_int[UMY]*s_int[UMY] + s_int[UMZ]*s_int[UMZ])/s_int[URHO])/s_int[URHO];
    // amrex::Real p_int = (parm.eos_gamma - 1.0)*s_int[URHO]*eint_int;
    // q_ext[QPRES] = p_int;
    // // T=Twall
    // amrex::Real T_int = p_int/(parm.Rspec*s_int[URHO]); 
    // q_ext[QT]    = max(prob_parm.Tw  +  dratio*(prob_parm.Tw - T_int),50.0);
    // // rho = eos(P,T)
    // q_ext[QRHO]  = q_ext[QPRES]/(parm.Rspec*q_ext[QT]);

    // // convert prims to cons
    // s_ext[URHO] = q_ext[QRHO];
    // s_ext[UMX] = q_ext[QRHO]*q_ext[QU];
    // s_ext[UMY] = q_ext[QRHO]*q_ext[QV];
    // s_ext[UMZ] = q_ext[QRHO]*q_ext[QW];
    // amrex::Real ekin_ext = 0.5*(q_ext[QU]*q_ext[QU] + q_ext[QV]*q_ext[QV] + q_ext[QW]*q_ext[QW]); 
    // amrex::Real eint_ext = q_ext[QPRES]/(q_ext[QRHO]*(parm.eos_gamma - 1.0));
    // s_ext[UET] = q_ext[QRHO]*(eint_ext + ekin_ext);
  }
}

} // namespace prob


#endif
