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
#include <Constants.h>

using namespace universal_constants;

namespace PROB {

//////////////////////////// Physical modelling ////////////////////////////////
struct ProbParm
{ 
  Real Ma    =  2.0_rt; 
  Real Re    =  10000.0_rt;       
  Real T_inf =  200.0_rt;           // K
  Real p_inf =  100000.0_rt;        // 1bar
  
  Real gam   = 1.4_rt;              // gamma
  Real Rgas  = gas_constant/(28.96_rt*mass_cgs2si);  // gas constant(air)

  // sound speed
  Real cs   = std::sqrt(gam*Rgas*T_inf);  
  
  // derived values  u,rho, energy
  Real u_inf = Ma*cs; 
  Real rho_inf = p_inf/(T_inf*Rgas);
  Real ekin    = rho_inf*u_inf*u_inf;
  Real eint    = p_inf/(gam -Real(1.0));

  // refinement region
  Real x0 = 1.0_rt;            
  Real y0 = 2.0_rt;
  Real z0 = 2.0_rt;

};

typedef closures_dt<indicies_t, visc_suth_t, cond_suth_t,
                    calorifically_perfect_gas_t<indicies_t>> ProbClosures;

typedef rhs_dt<rusanov_t<ProbClosures>, no_diffusive_t, no_source_t>  ProbRHS;

//typedef rhs_dt<skew_t<true,false, 4, ProbClosures>, no_diffusive_t, no_source_t > ProbRHS;


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
  Real P,  T, rho, ux;

  Real x = prob_lo[0] + (i+0.5_rt)*dx[0];
  Real y = prob_lo[1] + (j+0.5_rt)*dx[1];
  Real z = prob_lo[2] + (k+0.5_rt)*dx[2];

  Real u_small = pparm.u_inf*0.2;
  Real xstart = 0.45_rt;
  Real xend   = pparm.x0;
  Real dis = xstart - xend;
  
  // flow initialisation
  if (x < xstart) {
    ux = pparm.u_inf;
    T = pparm.T_inf;
    P = pparm.p_inf;
  }
  else if (pow(x-pparm.x0,2) + pow(y-pparm.y0,2) + pow(z-pparm.z0,2) < pow(0.6_rt,2) ) { 
    Real xx = x - xstart;
    ux = u_small;
    T = pparm.T_inf*3;
    P = pparm.p_inf;
  } else {
    ux = pparm.u_inf;
    T = pparm.T_inf*3;
    P = pparm.p_inf;
  }
  //--------------------------------

  rho = P/(closures.Rspec*T);


  // std::cout << "P = " << P << std::endl;
  // std::cout << "T = " << T << std::endl;
  // std::cout << "rho = " << rho << std::endl;

  // std::cout << "P = " << pparm.p_inf << std::endl;
  // std::cout << "T = " << pparm.T_inf << std::endl;
  // std::cout << "rho = " << pparm.rho_inf << std::endl;
  //  std::cout << "c = " << pparm.cs << std::endl;
  //  std::cout << "e = " << rho*closures.cv*T<< std::endl;
  //  std::cout << "eint = " << pparm.eint << std::endl;
   


 //exit(0);

  
  state(i,j,k,ProbClosures::URHO ) = rho;
  state(i,j,k,ProbClosures::UMX  ) = rho*ux;
  state(i,j,k,ProbClosures::UMY  ) = Real(0.0);
  state(i,j,k,ProbClosures::UMZ  ) = Real(0.0);
  state(i,j,k,ProbClosures::UET)   = rho*closures.cv*T+ 0.5_rt*rho*ux*ux;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE 
void user_tagging(int i, int j, int k, int nt, auto& tagfab, const auto &sdatafab, 
                  const Array4<bool>&ibfab, const auto& geomdata, 
                  const ProbParm& pparm , int level) {

  const Real* dx  = geomdata.CellSize();
  const Real x = (i+0.5_rt)*dx[0];
  const Real y = (j+0.5_rt)*dx[1];
  const Real z = (k+0.5_rt)*dx[2];
  Real xrel[3];
  // coordinate relative to object
  xrel[0]= x-pparm.x0; xrel[1]= y-pparm.y0; xrel[2]= z-pparm.z0;
  Real radius =xrel[0]*xrel[0] + xrel[1]*xrel[1] + xrel[2]*xrel[2];
  const Real Rmax = 0.5_rt*0.5_rt;const Real Rmin = 0.15_rt*0.15_rt;
  // initialize thresholds at all levels
  Real rhofluc_threshold[6] = {0.3_rt,0.6_rt,0.9_rt,1000_rt,1000_rt,1000_rt};

  // refinement first step
  if ( nt==0) {
    
    if (level==0 ) {    
      tagfab(i,j,k) = (radius < Rmax ) && (radius > Rmin);
    }

  }
  else
  {
    int URHO = ProbClosures::URHO; 
    // refine close to grads of density 
    Real drhox = std::abs(sdatafab(i+1,j,k,URHO) - sdatafab(i-1,j,k,URHO));
    Real drhoy = std::abs(sdatafab(i,j+1,k,URHO) - sdatafab(i,j-1,k,URHO));
    Real rhop  = sdatafab(i,j,k,URHO);
    Real rhofluc = std::sqrt(drhox*drhox + drhoy*drhoy)/rhop ;

    tagfab(i,j,k) = (rhofluc > rhofluc_threshold[level]);

    // always refine close to body (at all levels)
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
