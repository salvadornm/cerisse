#ifndef CNS_PROB_H_
#define CNS_PROB_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_ParmParse.H>
#include <Closures.h>
#include <RHS.h>


// 2D Shock Reflection

using namespace amrex;

namespace PROB {

// problem parameters 
struct ProbParm {
  Real gam = 1.4;
  Real p_l    = 116.5;    
  Real u_l    =  8.25  * cos(30.0 / 180.0 * M_PI);
  Real v_l    =  -8.25 * sin(30.0 / 180.0 * M_PI);
  Real rho_l  = 8.0; 
  Real eint_l =  p_l/ (gam - Real(1.0));
  Real c_l    = sqrt(gam*p_l/rho_l);

  Real p_r    = 1.0;    
  Real u_r    = 0.0;
  Real v_r    = 0.0;
  Real rho_r  = 1.4; 
  Real eint_r =   p_r/ (gam - Real(1.0));;
  Real c_r    = sqrt(gam*p_r/rho_r);

};

inline Vector<std::string> cons_vars_names={"Xmom","Ymom","Zmom","Energy","Density"};
inline Vector<int> cons_vars_type={1,2,3,0,0};

typedef closures_dt<indicies_t, visc_suth_t, cond_suth_t,
                    calorifically_perfect_gas_t<indicies_t>>
    ProbClosures;

template <typename cls_t > class user_source_t;


typedef rhs_dt<riemann_t<false, ProbClosures>, no_diffusive_t, no_source_t > ProbRHS;
//typedef rhs_dt<skew_t<false,false, 4, ProbClosures>, no_diffusive_t, no_source_t > ProbRHS;


void inline inputs() {
  ParmParse pp;

  pp.add("cns.order_rk", 3);   // -2, 1, 2 or 3"
  pp.add("cns.stages_rk", 3);  // 1, 2 or 3
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
  Real y = prob_lo[1] + (j + Real(0.5)) * dx[1];
  
  // local vars
  Real rhot,eint, T, Pt; 
  Real u[2];

  // final state

  if (x < 1. / 6. + 1. / sqrt(3.) * y) {
    rhot =  prob_parm.rho_l;
    u[0] =  prob_parm.u_l; u[1] =  prob_parm.v_l;
    Pt   =  prob_parm.p_l;    
    eint =  prob_parm.eint_l;
  } else {
    rhot =  prob_parm.rho_r;
    u[0] =  prob_parm.u_r; u[1] =  prob_parm.v_r;
    Pt    =  prob_parm.p_r;    
    eint =  prob_parm.eint_r;
  }

  //std::cout << " gam=" << cls.gamma << std::endl;
  // che
  eint =  Pt / (cls.gamma - Real(1.0));

  
  state(i, j, k, cls.URHO) = rhot;
  state(i, j, k, cls.UMX)  = rhot * u[0];
  state(i, j, k, cls.UMY)  = rhot * u[1];
  state(i, j, k, cls.UMZ)  = Real(0.0);  
  state(i, j, k, cls.UET)  = eint + Real(0.5) * rhot * (u[0] * u[0] + u[1] * u[1]);    
}

/////////////////////////////// BC /////////////////////////////////////////////
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
bcnormal(const Real x[AMREX_SPACEDIM], Real dratio, const Real s_int[ProbClosures::NCONS],
         const Real s_refl[ProbClosures::NCONS], Real s_ext[ProbClosures::NCONS], const int idir,
         const int sgn, const Real time, GeometryData const & /*geomdata*/,
         ProbClosures const &closures, ProbParm const &prob_parm) {

  const int URHO = ProbClosures::URHO;
  const int UMX  = ProbClosures::UMX;
  const int UMY  = ProbClosures::UMY;
  const int UMZ  = ProbClosures::UMZ;
  const int UET  = ProbClosures::UET;
   

  // x low
  if (idir == 0 && sgn == 1) { 
    // post shock conditions
    s_ext[URHO] = prob_parm.rho_l;
    s_ext[UMX]  = prob_parm.rho_l * prob_parm.u_l;
    s_ext[UMY]  = prob_parm.rho_l * prob_parm.v_l;
    s_ext[UMZ]  = 0.0;    
    s_ext[UET]  =   prob_parm.eint_l + 0.5 *prob_parm.rho_l*
      (prob_parm.u_l * prob_parm.u_l + prob_parm.v_l * prob_parm.v_l);
  }
  // x top
  if (idir == 0 && sgn == -1) { 
    // post shock conditions
    s_ext[URHO] = prob_parm.rho_l;
    s_ext[UMX]  = prob_parm.rho_l * prob_parm.u_l;
    s_ext[UMY]  = prob_parm.rho_l * prob_parm.v_l;
    s_ext[UMZ]  = 0.0;    
    s_ext[UET]  =   prob_parm.eint_l + 0.5 *prob_parm.rho_l*
      (prob_parm.u_l * prob_parm.u_l + prob_parm.v_l * prob_parm.v_l);  }
  // y low
  if (idir == 1 && sgn == 1) { 
    if (x[0] < 1. / 6.) {
      // post-shock conditions
      s_ext[URHO] = prob_parm.rho_l;
      s_ext[UMX]  = prob_parm.rho_l * prob_parm.u_l;
      s_ext[UMY]  = prob_parm.rho_l * prob_parm.v_l;
      s_ext[UMZ]  = 0.0;
    }
    else {
      // slip wall
      s_ext[URHO] = s_int[URHO];
      s_ext[UMX]  = s_int[UMX];
      s_ext[UMY]  = -s_int[UMY];
      s_ext[UMZ]  = s_int[UMZ];
      s_ext[UET]  = s_int[UET];
    }
  }
  // y top
  if (idir == 1 && sgn == -1) { 
    Real us = 10. * prob_parm.c_r / sqrt(3.) * 2.; // shock speed
    if (x[0] - us * time < 1. / 6. + 1. / sqrt(3.) * x[1]) {
      s_ext[URHO] = prob_parm.rho_l;
      s_ext[UMX] = prob_parm.rho_l * prob_parm.u_l;
      s_ext[UMY] = prob_parm.rho_l * prob_parm.v_l;
      s_ext[UMZ] = 0.0; 
      s_ext[UET]  =   prob_parm.eint_l + 0.5 *prob_parm.rho_l*
      (prob_parm.u_l * prob_parm.u_l + prob_parm.v_l * prob_parm.v_l);
    }
    else {
      s_ext[URHO] = prob_parm.rho_r;
      s_ext[UMX] = 0.0;
      s_ext[UMY] = 0.0;
      s_ext[UMZ] = 0.0;
      s_ext[UET] = prob_parm.eint_r;
    }
  
  } 
   
  //-----

}
///////////////////////////////AMR//////////////////////////////////////////////
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
user_tagging(int i, int j, int k, int nt_level, auto &tagfab,
             const auto &sdatafab, const auto &geomdata,
             const ProbParm &prob_parm, int level) {


  Real rhot = sdatafab(i,j,k,ProbClosures::URHO);

  Real dengrad_threshold = 0.5;
  Real drhox = Math::abs(sdatafab(i+1,j,k,ProbClosures::URHO) -
   sdatafab(i,j,k,ProbClosures::URHO))/rhot;

  Real drhoy = Math::abs(sdatafab(i,j+1,k,ProbClosures::URHO) -
   sdatafab(i,j-1,k,ProbClosures::URHO))/rhot;

  if (nt_level > 0)
  {
    // tag cells based on density values
    switch (level)
    {
      case 0:
        tagfab(i,j,k) = (rhot > 1.1 && rhot < 1.9);
        break;
      case 1:
        tagfab(i,j,k) = (rhot > 1.2 && rhot < 1.8);
        break;
      default:
        tagfab(i,j,k) = (rhot > 1.3 && rhot < 1.7);
        break;
    }
  }
}
////////////////////////////////////////////////////////////////////////////////




///////////////////////////////////////////////////////////////////////////////

} // namespace PROB
#endif
