#ifndef CNS_PROB_H_
#define CNS_PROB_H_

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_ParmParse.H>
#include <AMReX_AmrLevel.H>
#include <Closures.h>

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
#define NIMPS 3 // number of image points per ghost point

namespace PROB {

//////////////////////////// Physical modelling ////////////////////////////////
struct ProbParm
{ 
  Real p_inf = 167.9;
  Real T_inf = 54.8;
  Real u_inf = 1405.59;

  Real x0 = 0.1_rt;
  Real y0 = 0.2_rt;
  Real z0 = 0.25_rt;
};

typedef closures_derived_t<visc_suth_t, cond_suth_t, calorifically_perfect_gas_t> ProbClosures;

//////////////////////////// Numerical modelling ///////////////////////////////
void inline inputs() {
  ParmParse pp;

  // Numerical operators
  //-1 = N/A (Incase of periodic)
  // 0 = Interior           3 = Symmetry
  // 1 = Inflow             4 = SlipWall
  // 2 = Outflow            5 = NoSlipWall
  // 6 = user defined
  pp.addarr("cns.lo_bc", std::vector<int>{2,0,0});
  pp.addarr("cns.hi_bc", std::vector<int>{2,0,0});
  pp.add   ("cns.order_rk", 3); // -2, 1, 2 or 3"
  pp.add   ("cns.stages_rk", 3); // 1, 2 or 3
  pp.add   ("cns.rhs_euler", 1); // 0=false, 1=true
  pp.add   ("cns.rhs_visc", 0); // 0=false, 1=true
  pp.add   ("cns.rhs_source", 0); // 0=false, 1=true
  pp.add   ("cns.flux_euler", 1); // 0=riemann solver, 1=KEEP/AD, 2=WENO5
  pp.add   ("cns.order_keep", 4); // Order of accuracy=2, 4 or 6"
  pp.add   ("cns.art_diss",1); // 0=none, 1=artificial dissipation
  pp.add   ("cns.screen_output", 2); // 0=quiet, 1=verbose
  pp.add   ("cns.verbose", 1); // 0=quiet, 1=verbose

  // ibm
  pp.add   ("ib.move",0); // 0=false, 1=true
  pp.add   ("ib.plot_surf",0); // 0=false, 1=true

}

//////////////////////////// Initial conditions ////////////////////////////////
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void prob_initdata (int i, int j, int k, amrex::Array4<amrex::Real> const& state, amrex::GeometryData const& geomdata, ProbClosures const& closures, ProbParm const& pparm) {
  using amrex::Real;
  const Real* prob_lo = geomdata.ProbLo();
  const Real* prob_hi = geomdata.ProbHi();
  const Real* dx      = geomdata.CellSize();
  Real P, rho, T, ux;

  Real x = prob_lo[0] + (i+0.5_rt)*dx[0];
  Real y = prob_lo[1] + (j+0.5_rt)*dx[1];
  Real z = prob_lo[2] + (k+0.5_rt)*dx[2];
  // if (x < prob_hi[0]/8) {
  //     P = pparm.p_l;
  //     rho = pparm.rho_l;
  //     ux = pparm.u_l*0.0_rt;
  // } else {
  //     P = pparm.p_l;
  //     rho = pparm.rho_l;
  //     ux = pparm.u_l;
  // }

  // if (pow(x-pparm.x0,2) + pow(y-pparm.y0,2) + pow(z-pparm.z0,2) < pow(0.4_rt,2) ){
  //     P = pparm.p_inf;
  //     rho = pparm.rho_inf*0.7;
  //     ux = pparm.u_inf*0.9;
  // } else {
  //     P = pparm.p_inf;
  //     rho = pparm.rho_inf;
  //     ux = pparm.u_inf;
  // }

  Real u_small = pparm.u_inf;
  Real xstart = 0.005_rt;
  Real xend   = pparm.x0;
  Real dis = xstart - xend;
  Real grad = (pparm.u_inf - u_small)/(dis);
  if (x < xstart) {
    ux = pparm.u_inf;
    T = pparm.T_inf;
    P = pparm.p_inf;
  }
  else if (x < 0.05) { 
    Real xx = x - xstart;
    // ux = pparm.u_inf + grad*xx;
    ux = u_small*0.66;
    T = pparm.T_inf*3;
    P = pparm.p_inf;
  } else {
    // ux = u_small;
    ux = pparm.u_inf*0.3;
    T = pparm.T_inf*3;
    P = pparm.p_inf;
  }


  rho = P/(closures.Rspec*T);
  state(i,j,k,URHO ) = rho;
  state(i,j,k,UMX  ) = rho*ux;
  state(i,j,k,UMY  ) = Real(0.0);
  state(i,j,k,UMZ  ) = Real(0.0);
  Real eint = rho*closures.cv*T;
  state(i,j,k,UET) = eint + Real(0.5)*rho*ux*ux;
}

//
AMREX_GPU_DEVICE AMREX_FORCE_INLINE 
void user_tagging(int i, int j, int k, int nt, auto& tagfab, const auto &sdatafab, const auto &ibfab, const auto& geomdata, const ProbParm& pparm , int level) {
  const Real* prob_lo = geomdata.ProbLo();
  const Real* dx  = geomdata.CellSize();

  Real x = prob_lo[0] + (i+0.5_rt)*dx[0];
  Real y = prob_lo[1] + (j+0.5_rt)*dx[1];
  Real z = prob_lo[2] + (k+0.5_rt)*dx[2];

  if ( nt==0) {
    Real dengrad_threshold = 0.3;
    amrex::Real drhox = amrex::Math::abs(sdatafab(i+1,j,k,UMX) - sdatafab(i,j,k,UMX))/sdatafab(i,j,k,UMX);
    // if (drhox > dengrad_threshold) {
    //   for (int ii = -1; ii <= 1; ii++) {
    //     for (int jj = -1; jj <= 1; jj++) {
    //       for (int kk = -1; kk <= 1; kk++) {
    //         tagfab(i+ii,j+jj,k+kk) = true;
    //       }
    //     }
    //   }
    // }
    if (level==0 ) {
      if (pow(x-pparm.x0,2) + pow(y-pparm.y0,2) + pow(z-pparm.z0,2) < pow(0.08_rt,2) &&  (pow(x-pparm.x0,2) + pow(y-pparm.y0,2) + pow(z-pparm.z0,2) > pow(0.01_rt,2) )){
        tagfab(i,j,k) = true;
      }
    }


  }
  if (level==0 && nt>0) {
    Real dengrad_threshold = 0.3;
    amrex::Real drhox = amrex::Math::abs(sdatafab(i+1,j,k,URHO) - sdatafab(i,j,k,URHO))/sdatafab(i,j,k,URHO);
    // if (drhox > dengrad_threshold) {
    //   tagfab(i,j,k) = true;
    // }
    if (ibfab(i,j,k,1)) {
      for (int ii = -1; ii <= 1; ii++) {
        for (int jj = -1; jj <= 1; jj++) {
          for (int kk = -1; kk <= 1; kk++) {
            tagfab(i+ii,j+jj,k+kk) = true;
          }
        }
      }
    }
  }
  else if (level==1 && nt>0) {
    // Real dengrad_threshold = 0.1;
    // amrex::Real drhox = amrex::Math::abs(sdatafab(i+1,j,k,URHO) - sdatafab(i,j,k,URHO))/sdatafab(i,j,k,URHO);
    // if (sdatafab(i,j,k,URHO) > 600) {
    //   tagfab(i,j,k) = true;
    // }
    if (ibfab(i,j,k,1)) {
      for (int ii = -1; ii <= 1; ii++) {
        for (int jj = -1; jj <= 1; jj++) {
          for (int kk = -1; kk <= 1; kk++) {
            tagfab(i+ii,j+jj,k+kk) = true;
          }
        }
      }
    }
  }
}

//////////////////////////// Boundary conditions ///////////////////////////////
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
 * @param pparm ProbParm data as defined in pparm.H and initialised in
 * amrex_probinit.
 * @sa CnsFillExtDir
 * @sa CnsFillExtDir::operator()
 */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
bcnormal(const amrex::Real x[AMREX_SPACEDIM], amrex::Real dratio, const amrex::Real s_int[NCONS],
         const amrex::Real s_refl[NCONS], amrex::Real s_ext[NCONS],
         const int idir, const int sgn, const amrex::Real time,
         amrex::GeometryData const& /*geomdata*/,  ProbClosures const& closures, ProbParm const& pparm)
{
  // if (idir == 1) { // ylo or yhi

    amrex::Abort("bcnormal not coded");

}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void user_source(int i, int j, int k, const auto& state, const auto& rhs, const ProbParm& lprobparm, ProbClosures const& closures, auto const dx) {
};

}
#endif