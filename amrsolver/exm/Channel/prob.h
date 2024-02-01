#ifndef CNS_PROB_H_
#define CNS_PROB_H_

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
// #include <AMReX_PROB_AMR_F.H>
#include <AMReX_ParmParse.H>
#include <Closures.h>

using namespace amrex;


namespace PROB {

////////////////////////////////EQUATIONS///////////////////////////////////////
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

// Numerical operators
void inline inputs() {
  ParmParse pp;

  pp.addarr("cns.lo_bc", std::vector<int>{0,6,0});
  pp.addarr("cns.hi_bc", std::vector<int>{0,6,0});
  pp.add   ("cns.order_rk", 3); // -2, 1, 2 or 3"
  pp.add   ("cns.stages_rk", 3); // 1, 2 or 3
  pp.add   ("cns.rhs_euler", 1); // 0=false, 1=true
  pp.add   ("cns.rhs_visc", 1); // 0=false, 1=true
  pp.add   ("cns.rhs_source", 1); // 0=false, 1=true
  pp.add   ("cns.flux_euler", 1); // 0=riemann solver, 1=KEEP/AD, 2=WENO5
  pp.add   ("cns.order_keep", 4); // Order of accuracy=2, 4 or 6"
  pp.add   ("cns.art_diss", 0); // 0=none, 1=artificial dissipation
  pp.add   ("cns.screen_output", 50); //
  pp.add   ("cns.verbose", 1); // 0=quiet, 1=verbose
}

////////////////////////////////////////////////////////////////////////////////

///////////////////////////////CLOSURES/////////////////////////////////////////
typedef closures_derived_t<visc_suth_t, cond_suth_t, calorifically_perfect_gas_t> ProbClosures;
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
    Real h       = 0.006845_rt ;  // channel half height (m)
    Real Re_tau  = 220.0_rt;//Re_tau = rho_w u_tau h/mu_w
    Real Re_b    = 3.0e3;     //bulk Re = rho_b ub h/mu_w
    Real Ma_b    = 1.5_rt; //originally 1.5    //mach number bulk  = ub/aw  (aw is speed of sound at wall)
    Real u_tau = 35.0_rt; // (m/s)
    Real Tw    = 500.0_rt; // wall temperature (K)
    // Real Mt0   = 0.10;           // initial turbulent mach number

    Real Rgas = 287.0_rt;// specific gas constant air (J/kgK)
    // rhow, Pw and P0
    Real muw = 2.670533348e-5; // mu at 500K for air
    Real rho_w = Re_tau * muw / (u_tau*h);
    Real Pw   = rho_w*Rgas*Tw; // csal_press(rho_w,T_w,nm_csal)

    // bulk velocity and density
    Real aw = pow(1.4_rt*Rgas*Tw,0.5_rt);
    Real ub = Ma_b*aw;
    Real rho_b = Re_b*muw/(ub*h);
    Real Tb = (Pw/Rgas*rho_b);

    // body force in x
    Real fx =  rho_w*u_tau*u_tau/(h*rho_b);
};

AMREX_GPU_DEVICE inline
void prob_initdata (int i, int j, int k, Array4<Real> const& state, GeometryData const& geomdata, ProbClosures const& cls, const ProbParm& pparm) {
  const Real* prob_lo = geomdata.ProbLo();
  const Real* prob_hi = geomdata.ProbHi();
  const Real* dx      = geomdata.CellSize();

  Real pi = 3.14159265358979323846_rt;
  Real x = prob_lo[0] + (i+Real(0.5))*dx[0];
  Real y = prob_lo[1] + (j+Real(0.5))*dx[1];
  Real z = prob_lo[2] + (k+Real(0.5))*dx[2];

  Real lx = prob_hi[0] - prob_lo[0];
  Real ly = prob_hi[1] - prob_lo[1];
  Real lz = prob_hi[2] - prob_lo[2];

  if (y > pparm.h) {y = 2*pparm.h-y;}

  // Laminar velocity and constant temperature
  // Real ux = 1.5*pparm.ub*(1 - (1/std::pow(pparm.h,2))*std::pow((y-pparm.h),2));
  // Real ux = 0.0_rt;

  // Real T = pparm.Tw+0.5*pparm.Tw*(1 - (1/std::pow(pparm.h,2))*std::pow((y-pparm.h),2));
  // Real T = pparm.Tw;

  // Turbulent profiles
  Real T = 677.04110513 + 255.46817*pow(y,0.5) - 176.38967315*exp(-3081.95957181*y)*1.025;
  Real ux = 581.54524091 + 2687.84321585*pow(y,0.5) - 593.73493364* exp(-2008.87100276*y);

  int nx = 1;
  int ny = 2;
  int nz = 3;

  // srand(10);
  // Print() << (Real)rand()/ RAND_MAX << std::endl;
  // amrex::RandomEngine re;
// Real RandomNormal (Real mean, Real stddev, RandomEngine const& random_engine)
  // Real rand_num = 1.0;//amrex::RandomNormal(0.0, 0.0,re);
#if AMREX_USE_GPU
  Real rand_num = 0.5;
#else
  Real rand_num = (Real)rand()/ RAND_MAX;
#endif

  Real A       = 100*2.5;
  Real u_prime = rand_num*A*cos(2*pi*x*nx/lx)*sin(pi*y*ny/ly)*cos(2*pi*z*nz/lz);
  Real v_prime = rand_num*-A*sin(2*pi*x*nx/lx)*cos(pi*y*ny/ly)*sin(2*pi*z*nz/lz);
  Real w_prime = rand_num*-A*sin(2*pi*x*nx/lx)*cos(pi*y*ny/ly)*sin(2*pi*z*nz/lz);

  Real rho =  pparm.Pw/(cls.Rspec*T);
  state(i,j,k,URHO ) = rho;
  state(i,j,k,UMX  ) = rho*(ux+u_prime);
  state(i,j,k,UMY  ) = rho*v_prime;
  state(i,j,k,UMZ  ) = rho*w_prime;
  state(i,j,k,UET  ) = rho*cls.cv*T + Real(0.5)*rho*(ux*ux + u_prime*u_prime + v_prime*v_prime + w_prime*w_prime);
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE 
void bcnormal(const Real x[AMREX_SPACEDIM], Real dratio, const Real s_int[NCONS], const Real s_refl[NCONS], Real s_ext[NCONS], const int idir, const int sgn, const Real time, GeometryData const& /*geomdata*/,  ProbClosures const& cls, ProbParm const& pparm)
{
  if (idir == 1) { // ylo or yhi

    Real q_ext[NPRIM] = {0.0};
    // no-slip
    q_ext[QU]    = -s_int[UMX]/s_int[URHO];
    q_ext[QV]    = -s_int[UMY]/s_int[URHO];
    q_ext[QW]    = -s_int[UMZ]/s_int[URHO];

    // dp/dn = 0
    amrex::Real eint_int = (s_int[UET] - 0.5*(s_int[UMX]*s_int[UMX] + s_int[UMY]*s_int[UMY] + s_int[UMZ]*s_int[UMZ])/s_int[URHO])/s_int[URHO];
    amrex::Real p_int = (cls.gamma - 1.0)*s_int[URHO]*eint_int;
    q_ext[QPRES] = p_int;
    // T=Twall
    amrex::Real T_int = p_int/(cls.Rspec*s_int[URHO]); 
    q_ext[QT]    = max(pparm.Tw  +  dratio*(pparm.Tw - T_int),50.0);
    // rho = eos(P,T)
    q_ext[QRHO]  = q_ext[QPRES]/(cls.Rspec*q_ext[QT]);

    // convert prims to cons
    s_ext[URHO] = q_ext[QRHO];
    s_ext[UMX] = q_ext[QRHO]*q_ext[QU];
    s_ext[UMY] = q_ext[QRHO]*q_ext[QV];
    s_ext[UMZ] = q_ext[QRHO]*q_ext[QW];
    amrex::Real ekin_ext = 0.5*(q_ext[QU]*q_ext[QU] + q_ext[QV]*q_ext[QV] + q_ext[QW]*q_ext[QW]); 
    amrex::Real eint_ext = q_ext[QPRES]/(q_ext[QRHO]*(cls.gamma - 1.0));
    s_ext[UET] = q_ext[QRHO]*(eint_ext + ekin_ext);
  }
}

AMREX_GPU_DEVICE inline
void user_source(int i, int j, int k, const auto& state, const auto& rhs, const ProbParm& pparm, ProbClosures const& cls, auto const dx) {
  rhs(i,j,k,UMX) +=  state(i,j,k,URHO)*pparm.fx;
  rhs(i,j,k,UET) +=  state(i,j,k,UMX )*pparm.fx;
}
////////////////////////////////////////////////////////////////////////////////


///////////////////////////////AMR//////////////////////////////////////////////
AMREX_GPU_DEVICE inline 
void user_tagging(int i, int j, int k, int nt_lev , auto& tagfab, const auto &sdatafab, const auto& geomdata, const ProbParm& pparm , int level) {

    int jmax = geomdata.Domain().bigEnd(1);
    int jmin = geomdata.Domain().smallEnd(1);

    const Real* prob_hi = geomdata.ProbHi();
    const Real* dx      = geomdata.CellSize();
    Real y = (j+Real(0.5))*dx[1];
    Real factor = 0.5;
    tagfab(i,j,k) = y < pparm.h*factor || y > prob_hi[1] - pparm.h*factor;
  }
////////////////////////////////////////////////////////////////////////////////
} // namespace PROB
#endif
