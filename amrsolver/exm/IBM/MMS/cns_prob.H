#ifndef CNS_PROB_H
#define CNS_PROB_H

#include <AMReX_ParmParse.H>
#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>
#include <Closures.h>

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

#define NGHOST 2

namespace PROB {

struct ProbParm
{
  const Real pi = 3.1415926535897932384626433_rt;
  Real Lambda = 1.0_rt; // length of the domain, assuming it is a cube
  const Real ampl_p = 1.0_rt;
  const Real ampl_t = 1.0_rt;
  const Real ampl_u = 1.0_rt;
  const Real ampl_v = 1.0_rt;
  const Real ampl_w = 1.0_rt;
  const Real n_p = 1.0_rt;
  const Real n_t = 1.0_rt;
  const Real n_u = 1.0_rt;
  const Real n_v = 1.0_rt;
  const Real n_w = 1.0_rt;

  const Real alpha_p = 2.0_rt*n_p*pi/Lambda;
  const Real alpha_t = 2.0_rt*n_t*pi/Lambda;
  const Real alpha_u = 2.0_rt*n_u*pi/Lambda;
  const Real alpha_v = 2.0_rt*n_v*pi/Lambda;
  const Real alpha_w = 2.0_rt*n_w*pi/Lambda;

  Real fu (Real x, Real y, Real z) const {
    return ampl_u * cos(alpha_u*x)*sin(alpha_u*y)*sin(alpha_u*z);
  }

  Real fv (Real x, Real y, Real z) const {
    return ampl_v * cos(alpha_v*x)*sin(alpha_v*y)*sin(alpha_v*z);
  }

  Real fw (Real x, Real y, Real z) const {
    return ampl_w * cos(alpha_w*x)*sin(alpha_w*y)*sin(alpha_w*z);
  }

  Real ft (Real x, Real y, Real z) const {
    return 4.0_rt + ampl_t * cos(alpha_t*x)*sin(alpha_t*y)*cos(alpha_t*z);
  }

  Real fp (Real x, Real y, Real z) const {
    return 4.0_rt + ampl_p * cos(alpha_p*x)*sin(alpha_p*y)*cos(alpha_p*z);
  }

  Real dudx (Real x, Real y, Real z) const {
    return -ampl_u * alpha_u * sin(alpha_u*x) * sin(alpha_u*y) * sin(alpha_u*z) ;}

  Real dudy (Real x, Real y, Real z) const {
    return ampl_u * alpha_u * cos(alpha_u*x) * cos(alpha_u*y) * sin(alpha_u*z) ;}
  
  Real dudz (Real x, Real y, Real z) const {
    return ampl_u * alpha_u * cos(alpha_u*x) * sin(alpha_u*y) * cos(alpha_u*z) ;}

  Real dvdx (Real x, Real y, Real z) const {
    return - ampl_v * alpha_v * sin(alpha_v*x) * sin(alpha_v*y) * sin(alpha_v*z) ;}

  Real dvdy (Real x, Real y, Real z) const {
    return  ampl_v * alpha_v * cos(alpha_v*x) * cos(alpha_v*y) * sin(alpha_v*z) ;}

  Real dvdz (Real x, Real y, Real z) const {
    return  ampl_v * alpha_v * cos(alpha_v*x) * sin(alpha_v*y) * cos(alpha_v*z) ;}

  Real dwdx (Real x, Real y, Real z) const {
    return - ampl_w * alpha_w * sin(alpha_w*x) * sin(alpha_w*y) * sin(alpha_w*z) ;}

  Real dwdy (Real x, Real y, Real z) const {
    return  ampl_w * alpha_w * cos(alpha_w*x) * cos(alpha_w*y) * sin(alpha_w*z) ;}

  Real dwdz (Real x, Real y, Real z) const {
    return  ampl_w * alpha_w * cos(alpha_w*x) * sin(alpha_w*y) * cos(alpha_w*z) ;}

  Real dpdx (Real x, Real y, Real z) const {
    return - ampl_p * alpha_p * sin(alpha_p*x) * sin(alpha_p*y) * cos(alpha_p*z) ;}

  Real dpdy (Real x, Real y, Real z) const {
    return ampl_p * alpha_p * cos(alpha_p*x) * cos(alpha_p*y) * cos(alpha_p*z) ;}

  Real dpdz (Real x, Real y, Real z) const {
    return - ampl_p * alpha_p * cos(alpha_p*x) * sin(alpha_p*y) * sin(alpha_p*z);}

  Real dTdx (Real x, Real y, Real z) const {
    return - ampl_t * alpha_t * sin(alpha_t*x) * sin(alpha_t*y) * cos(alpha_t*z) ;}

  Real dTdy (Real x, Real y, Real z) const {
    return ampl_t * alpha_t * cos(alpha_t*x) * cos(alpha_t*y) * cos(alpha_t*z) ;}

  Real dTdz (Real x, Real y, Real z) const {
    return - ampl_t * alpha_t * cos(alpha_t*x) * sin(alpha_t*y) * sin(alpha_t*z);}

};

///////////////////////////////CLOSURES/////////////////////////////////////////
typedef closures_derived_base_t<visc_suth_t, cond_suth_t, calorifically_perfect_gas_t, ProbParm> ProbClosures;
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
  // pp.addarr("cns.lo_bc", std::vector<int>{2,-1,-1});
  // pp.addarr("cns.hi_bc", std::vector<int>{2,-1,-1});
  // pp.add   ("cns.order_rk", 2); // -2, 1, 2 or 3"
  // pp.add   ("cns.stages_rk", 2); // 1, 2 or 3
  // pp.add   ("cns.rhs_euler", 1); // 0=false, 1=true
  // pp.add   ("cns.rhs_visc", 0); // 0=false, 1=true
  // pp.add   ("cns.rhs_source", 0); // 0=false, 1=true
  // pp.add   ("cns.flux_euler", 0); // 0=riemann solver, 1=KEEP/AD, 2=WENO5
  // pp.add   ("cns.order_keep", 4); // Order of accuracy=2, 4 or 6"
  // pp.add   ("cns.art_diss", 0); // 0=none, 1=artificial dissipation
  // pp.add   ("cns.nghost",2);
  // pp.add   ("cns.screen_output", 1); // 0=quiet, 1=verbose
  // pp.add   ("cns.verbose", 1); // 0=quiet, 1=verbose

  // lets add viscosity type and different flux dissipations!


  // debugging
  // pp.add("amrex.fpe_trap_invalid",1);
  // pp.add("amrex.fpe_trap_zero",1);
  // pp.add("amrex.fpe_trap_overflow",1);

}



AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void prob_initdata (int i, int j, int k, amrex::Array4<amrex::Real> const& state, amrex::GeometryData const& geomdata, ProbClosures const& cls, ProbParm const& pparm) {
  // Geometry
  const Real* prob_lo = geomdata.ProbLo();
  const Real* prob_hi = geomdata.ProbHi();
  const Real* dx      = geomdata.CellSize();
  const Real x = prob_lo[0] + (i + 0.5) * dx[0];
  const Real y = prob_lo[1] + (j + 0.5) * dx[1];
  const Real z = prob_lo[2] + (k + 0.5) * dx[2];

  // MMS functions
  Real u[3],p,T,rho,eint;
  u[0] = pparm.fu(x,y,z);
  u[1] = pparm.fv(x,y,z);
  u[2] = pparm.fw(x,y,z);
  p = pparm.fp(x,y,z);
  T = pparm.ft(x,y,z);
  rho  = p/(cls.Rspec*T); 
  eint = cls.cv*T;

  // Set the state
  state(i, j, k, URHO) = rho;
  state(i, j, k, UMX)  = rho * u[0];
  state(i, j, k, UMY)  = rho * u[1];
  state(i, j, k, UMZ)  = rho * u[2];
  state(i, j, k, UET)  = rho * (eint + Real(0.5) * (u[0]*u[0] + u[1]*u[1] + u[2]*u[2]));
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE 
void user_tagging(int i, int j, int k, auto& tagfab, const auto &sdatafab, const auto &ibfab, const auto& geomdata, const ProbParm& prob_parm , int level) {

      Real dengrad_threshold = 0.5;
      amrex::Real drhox = amrex::Math::abs(sdatafab(i+1,j,k,URHO) - sdatafab(i-1,j,k,URHO))/sdatafab(i,j,k,URHO);
      if (drhox > dengrad_threshold) {
        tagfab(i,j,k) = true;
        tagfab(i+1,j,k) = true;
        tagfab(i+2,j,k) = true;
        tagfab(i+3,j,k) = true;
      }
  }

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
bcnormal(const amrex::Real x[AMREX_SPACEDIM], amrex::Real dratio, const amrex::Real s_int[NCONS],
         const amrex::Real s_refl[NCONS], amrex::Real s_ext[NCONS],
         const int idir, const int sgn, const amrex::Real time,
         amrex::GeometryData const& /*geomdata*/,  ProbClosures const& closures, ProbParm const& prob_parm)
{
  amrex::Abort("bcnormal not defined");
}


AMREX_GPU_DEVICE inline
void user_source(int i, int j, int k, const auto& state, const auto& rhs, const ProbParm& pparm, ProbClosures const& cls, auto const dx) {


  const Real x = (i + 0.5_rt) * dx[0];
  const Real y = (j + 0.5_rt) * dx[1];
  const Real z = (k + 0.5_rt) * dx[2];

  Real rho = 1/cls.Rspec;

  // rhs(i,j,k,URHO) =rhs(i,j,k,URHO) + rho*(pparm.dudx(x,y,z) + pparm.dvdy(x,y,z) + pparm.dwdz(x,y,z));

  Real error = rhs(i,j,k,URHO) + rho*(pparm.dudx(x,y,z) + pparm.dvdy(x,y,z) + pparm.dwdz(x,y,z));
  rhs(i,j,k,URHO) = error;

  // rho (duu/dx + duv/dy + duw/dz) assuming rho=constant
  error = rhs(i,j,k,UMX) + rho*( 2*pparm.fu(x,y,z)*pparm.dudx(x,y,z) 
                 + pparm.fv(x,y,z)*pparm.dudy(x,y,z) + pparm.fu(x,y,z)*pparm.dvdy(x,y,z)
                 + pparm.fw(x,y,z)*pparm.dudz(x,y,z) + pparm.fu(x,y,z)*pparm.dwdz(x,y,z)) 
                 + pparm.dpdx(x,y,z);
  rhs(i,j,k,UMX) = error;

  // rho (dvu/dx + dvv/dy + dvw/dz) assuming rho=constant
  error = rhs(i,j,k,UMY) + rho*( pparm.fu(x,y,z)*pparm.dvdx(x,y,z) + pparm.fv(x,y,z)*pparm.dudx(x,y,z) 
                 + 2*pparm.fv(x,y,z)*pparm.dvdy(x,y,z) 
                 + pparm.fw(x,y,z)*pparm.dvdz(x,y,z) + pparm.fv(x,y,z)*pparm.dwdz(x,y,z)) 
                 + pparm.dpdy(x,y,z);
  rhs(i,j,k,UMY) = error;

  // rho (dwu/dx + dwv/dy + dww/dz) assuming rho=constant
  error = rhs(i,j,k,UMZ) + rho*(pparm.fu(x,y,z)*pparm.dwdx(x,y,z) + pparm.fw(x,y,z)*pparm.dudx(x,y,z) 
        + pparm.fv(x,y,z)*pparm.dwdy(x,y,z) + pparm.fw(x,y,z)*pparm.dvdy(x,y,z)
        + 2*pparm.fw(x,y,z)*pparm.dwdz(x,y,z))
        + pparm.dpdz(x,y,z);
  rhs(i,j,k,UMZ) = error;

  // d (rho u ht) = d(rho u (et + p/rho)) = d(rho u et + Pu) = rho d(u et + Pu/rho) = rho [ d(u et) + 1/rho d(Pu) ] = rho [ et d(u) + u d(et) + (1/rho) ( u d(P) + P d(u) )]

  // et = e + 1/2 (u^2 + v^2 + w^2) = cv T + 1/2 (u^2 + v^2 + w^2)
  Real et = cls.cv*pparm.ft(x,y,z) + 0.5_rt*(pparm.fu(x,y,z)*pparm.fu(x,y,z) + pparm.fv(x,y,z)*pparm.fv(x,y,z) + pparm.fw(x,y,z)*pparm.fw(x,y,z));

  // x-direction: rho [ et d(u)/dx + u d(et)/dx + (1/rho) ( u d(P)dx + P d(u)dx )]
  // d(et)/dx = cv dT/dx + 2*0.5(udu/dx + vdv/dx + wdw/dx)
  Real det = cls.cv*pparm.dTdx(x,y,z) + (pparm.fu(x,y,z)*pparm.dudx(x,y,z) + pparm.fv(x,y,z)*pparm.dvdx(x,y,z) + pparm.fw(x,y,z)*pparm.dwdx(x,y,z));

  error = rho*( et*pparm.dudx(x,y,z) +  pparm.fu(x,y,z)*det + (1/rho)*(pparm.fu(x,y,z)*pparm.dpdx(x,y,z) + pparm.fp(x,y,z)*pparm.dudx(x,y,z)) );

  // y-direction: rho [ et d(v) + v d(et) + (1/rho) ( v d(P) + P d(v) )]
  det = cls.cv*pparm.dTdy(x,y,z) + (pparm.fu(x,y,z)*pparm.dudy(x,y,z) + pparm.fv(x,y,z)*pparm.dvdy(x,y,z) + pparm.fw(x,y,z)*pparm.dwdy(x,y,z));

  error += rho*( et*pparm.dvdy(x,y,z) +  pparm.fv(x,y,z)*det + (1/rho)*(pparm.fv(x,y,z)*pparm.dpdy(x,y,z) + pparm.fp(x,y,z)*pparm.dvdy(x,y,z)) );

  // z-direction: rho [ et d(w) + w d(et) + (1/rho) ( w d(P) + P d(w) )]
  det = cls.cv*pparm.dTdz(x,y,z) + (pparm.fu(x,y,z)*pparm.dudz(x,y,z) + pparm.fv(x,y,z)*pparm.dvdz(x,y,z) + pparm.fw(x,y,z)*pparm.dwdz(x,y,z));

  error += rho*( et*pparm.dwdz(x,y,z) +  pparm.fw(x,y,z)*det + (1/rho)*(pparm.fw(x,y,z)*pparm.dpdz(x,y,z) + pparm.fp(x,y,z)*pparm.dwdz(x,y,z)) );

  // Print() << x << " " << y << " " << z << " " << rhs(i,j,k,UET)/2 << " " << error << std::endl;
  // exit(0);

  error = error + rhs(i,j,k,UET);

  rhs(i,j,k,UET) = error;

}

}

#endif
