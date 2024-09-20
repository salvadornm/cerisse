#ifndef CNS_PROB_H_
#define CNS_PROB_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_ParmParse.H>
#include <Closures.h>
#include <RHS.h>


using namespace amrex;

namespace PROB {

// problem parameters  (1:bottom   2:top)
struct ProbParm {
  Real p_int = 2.0;
  Real rho_1 = 1.0;
  Real rho_2 = 2.0;
  Real grav = -1.0; 
  Real eps =  0.025;
};

inline Vector<std::string> cons_vars_names={"Xmom","Ymom","Zmom","Energy","Density"};
inline Vector<int> cons_vars_type={1,2,3,0,0};

typedef closures_dt<indicies_t, visc_suth_t, cond_suth_t,
                    calorifically_perfect_gas_t<indicies_t>>
    ProbClosures;

template <typename cls_t > class user_source_t;



typedef rhs_dt<riemann_t<false, ProbClosures>, no_diffusive_t, user_source_t<ProbClosures>>
    ProbRHS;


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
  
  Real Pt, rhot, uxt,uyt;
  Real Lint = prob_hi[1] / 2; // half-domain lenght in y
  Real Pint = prob_parm.p_int; // interface Pressure

  const Real freq = Real(8)*Real(3.14159265359); // wavelength = x-domain

  Real yrel = y - Lint;
  Real delta  = dx[1]/5;  // transition region between top/bottom (if < dx, sharp)
  Real step = Real(0.5) + Real(0.5)*tanh(yrel/delta);

  rhot = step*prob_parm.rho_2 + (Real(1.0) -step)*prob_parm.rho_1;

  Pt = Pint + rhot*prob_parm.grav*yrel; // hydrostatic pressure
  
  Real csound = sqrt(cls.gamma*Pt/rhot); 

  uxt = Real(0.0);
  uyt = -prob_parm.eps*cos(freq*x)*csound; // perturbation in y-component
  
  state(i, j, k, cls.URHO) = rhot;
  state(i, j, k, cls.UMX)  = rhot * uxt;
  state(i, j, k, cls.UMY)  = rhot * uyt;
  state(i, j, k, cls.UMZ)  = Real(0.0);

  Real et = Pt / (cls.gamma - Real(1.0));
  state(i, j, k, cls.UET) = et + Real(0.5) * rhot * (uxt * uxt + uyt * uyt); 

}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
bcnormal(const Real x[AMREX_SPACEDIM], Real dratio, const Real s_int[5],
         const Real s_refl[ProbClosures::NCONS], Real s_ext[5], const int idir,
         const int sgn, const Real time, GeometryData const & /*geomdata*/,
         ProbClosures const &closures, ProbParm const &prob_parm) {
}

///////////////////////////////SOURCE TERM /////////////////////////////////////
template <typename cls_t>
class user_source_t {
  public:
  void inline src(const amrex::MFIter &mfi,
                  const amrex::Array4<const amrex::Real> &prims,
                  const amrex::Array4<amrex::Real> &rhs, const cls_t *cls_d,
                  amrex::Real dt){

    const Box bx = mfi.tilebox();

    ProbParm const prob_parm; 

    amrex::ParallelFor(bx,
      [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      { 
      const auto& cls = *cls_d;
    
    //  user_source(i,j,k,prims,rhs,lprobparm,cls,dx);           
      Real rho = prims(i, j, k, cls.QRHO);
      rhs(i,j,k,cls.UMY) += prob_parm.grav*rho; 
      rhs(i,j,k,cls.UET) += prob_parm.grav*rho*prims(i, j, k, cls.QV);
      });

  };
};

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

} // namespace PROB
#endif
