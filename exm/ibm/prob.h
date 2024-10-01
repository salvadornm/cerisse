#ifndef CNS_PROB_H_
#define CNS_PROB_H_

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_ParmParse.H>
#include <AMReX_AmrLevel.H>
#include <Closures.h>
#include <RHS.h>
#include <eib.h>
#include <ratio>

namespace PROB {

//////////////////////////// Physical modelling ////////////////////////////////
struct ProbParm
{ 
  Real p_inf = 1.0e5;
  Real T_inf = 300.0;
  Real u_inf = 200.0; // 1200

  Real x0 = 2.0_rt;  
  Real y0 = 2.0_rt;
  Real z0 = 2.0_rt;

  // derived values
  Real gam = 1.4;
  Real eint    = p_inf/(gam -Real(1.0));
  Real rho_inf = p_inf*29.0_rt/(T_inf*8.314_rt*1000_rt);
  Real ekin = rho_inf*u_inf*u_inf;

};


typedef closures_dt<indicies_t, visc_suth_t, cond_suth_t,
                    calorifically_perfect_gas_t<indicies_t>>
    ProbClosures;

typedef rhs_dt<rusanov_t<ProbClosures>, no_diffusive_t, no_source_t>  ProbRHS;

//typedef rhs_dt<skew_t<false,false, 2, ProbClosures>, no_diffusive_t, no_source_t > ProbRHS;


typedef std::ratio<5,5> d_image;
typedef eib_t<1,1,d_image,ProbClosures> ProbIB;

void inline inputs() {
  ParmParse pp;

  pp.add("cns.order_rk", 3);   // -2, 1, 2 or 3"
  pp.add("cns.stages_rk", 3);  // 1, 2 or 3

  pp.add   ("ib.move",0); // 0=false, 1=true
  pp.add   ("ib.plot_surf",1); // 0=false, 1=true
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

  Real u_small = pparm.u_inf*0.2;
  Real xstart = 0.45_rt;
  Real xend   = pparm.x0;
  Real dis = xstart - xend;
  Real grad = (pparm.u_inf - u_small)/(dis);
  
  // if (x < xstart) {
  //   ux = pparm.u_inf;
  //   T = pparm.T_inf;
  //   P = pparm.p_inf;
  // }
  // else if (pow(x-pparm.x0,2) + pow(y-pparm.y0,2) + pow(z-pparm.z0,2) < pow(0.6_rt,2) ) { 
  //   Real xx = x - xstart;
  //   // ux = pparm.u_inf + grad*xx;
  //   ux = u_small;
  //   T = pparm.T_inf*3;
  //   P = pparm.p_inf;
  // } else {
  //   ux = pparm.u_inf;
  //   T = pparm.T_inf*3;
  //   P = pparm.p_inf;
  // }

  ux = Real(0.0);//pparm.u_inf;
  T = pparm.T_inf;
  P = pparm.p_inf;
  //--------------------------------


  rho = P/(closures.Rspec*T);
  state(i,j,k,ProbClosures::URHO ) = rho;
  state(i,j,k,ProbClosures::UMX  ) = rho*ux;
  state(i,j,k,ProbClosures::UMY  ) = Real(0.0);
  state(i,j,k,ProbClosures::UMZ  ) = Real(0.0);
  Real eint = rho*closures.cv*T;
  state(i,j,k,ProbClosures::UET) = eint + Real(0.5)*rho*ux*ux;
}

//
AMREX_GPU_DEVICE AMREX_FORCE_INLINE 
void user_tagging(int i, int j, int k, int nt, auto& tagfab, const auto &sdatafab, const Array4<bool>&ibfab, const auto& geomdata, const ProbParm& pparm , int level) {

  const Real* dx  = geomdata.CellSize();
  Real x = (i+0.5_rt)*dx[0];
  Real y = (j+0.5_rt)*dx[1];
  Real z = (k+0.5_rt)*dx[2];

  if ( nt==0) {
    Real dengrad_threshold = 0.3;
    amrex::Real drhox = amrex::Math::abs(sdatafab(i+1,j,k,ProbClosures::UMX) - sdatafab(i,j,k,ProbClosures::UMX))/sdatafab(i,j,k,ProbClosures::UMX);
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
      if (pow(x-pparm.x0,2) + pow(y-pparm.y0,2) + pow(z-pparm.z0,2) < pow(0.50_rt,2) &&  (pow(x-pparm.x0,2) + pow(y-pparm.y0,2) + pow(z-pparm.z0,2) > pow(0.15_rt,2) )){
        tagfab(i,j,k) = true;
      }
    }


  }
  if (level==0 && nt>0) {
    Real dengrad_threshold = 0.3;
    amrex::Real drhox = amrex::Math::abs(sdatafab(i+1,j,k,ProbClosures::URHO) - sdatafab(i,j,k,ProbClosures::URHO))/sdatafab(i,j,k,ProbClosures::URHO);
    if (drhox > dengrad_threshold) {
      tagfab(i,j,k) = true;
    }
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
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
bcnormal(const amrex::Real x[AMREX_SPACEDIM], amrex::Real dratio, const amrex::Real s_int[ProbClosures::NCONS],
         const amrex::Real s_refl[ProbClosures::NCONS], amrex::Real s_ext[ProbClosures::NCONS],
         const int idir, const int sgn, const amrex::Real time,
         amrex::GeometryData const& /*geomdata*/,  ProbClosures const& closures, ProbParm const& pparm)  
{
  const int URHO = ProbClosures::URHO;
  const int UMX  = ProbClosures::UMX;
  const int UMY  = ProbClosures::UMY;
  const int UMZ  = ProbClosures::UMZ;
  const int UET  = ProbClosures::UET;
  const int face = (idir+1)*sgn;

  switch(face)
  {
    case  -1:  // EAST
      
      s_ext[URHO] = pparm.rho_inf;
      s_ext[UMX]  = pparm.rho_inf*pparm.u_inf;
      s_ext[UMY]  = Real(0.0);    
      s_ext[UMZ]  = Real(0.0);
      s_ext[UET]  = pparm.eint + pparm.ekin;
  
      break;
    case   1:  //WEST
      break;      
    default:

      break;
  }


}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void user_source(int i, int j, int k, const Array4<Real>& state, const auto& rhs, const ProbParm& lprobparm, ProbClosures const& closures, auto const dx) {
};

}
#endif
