#ifndef walltypes_H_
#define walltypes_H_

#include <AMReX_EB2.H>
#include <AMReX_EB2_IF.H>
#include <AMReX_ParmParse.H>

#include <algorithm>
#include <cmath>

#include <EBMultiFab.h>

#ifdef USE_PELEPHYSICS
#include <TransPele.h>
#endif

////////////////////////////////////////////////////////////////////////////
// test wall
template <typename cls_t>
class test_wall_t
{
  public:

    test_wall_t() {}

    ~test_wall_t() {}

    static void inline wall_flux(const auto &geomdata, int i, int j, int k, const Real norm[AMREX_SPACEDIM], 
      amrex::GpuArray<amrex::Real, cls_t::NPRIM>& prims, amrex::GpuArray<amrex::Real, cls_t::NCONS>& fluxw,const cls_t* cls) {

        printf("TEST wall test passed\n ");
        // only set-non zero values (by default fluxw=0)
                
        // new tangential vectors 
        Real tan1[AMREX_SPACEDIM];
        tan1[0]  =  norm[1];tan1[1]  = -norm[0];
#if AMREX_SPACEDIM==3        
        tan1[2] = 0.0;
        Real tan2[AMREX_SPACEDIM];
        // tan2 = cross_productv(tan1,norm)
        tan2[0] = tan1[1] * norm[2] - tan1[2] * norm[1]; // x component
        tan2[1] = tan1[2] * norm[0] - tan1[0] * norm[2]; // y component
        tan2[2] = tan1[0] * norm[1] - tan1[1] * norm[0];  // z component      
#endif        
        // change of coordinates x,y,z  --> n,t1,t2
        Real vel[AMREX_SPACEDIM]; Real velxyz[AMREX_SPACEDIM];        
        vel[0] = 0.0;
        vel[1] = prims[cls_t::QU]*tan1[0] + prims[cls_t::QV]*tan1[1];
#if AMREX_SPACEDIM==3  
        vel[2] = prims[cls_t::QU]*tan2[0] + prims[cls_t::QV]*tan2[1] + prims[cls_t::QW]*tan2[2] ;                
#endif

        // change of coordinates  n,t1,t2 --> x,y,z 
#if AMREX_SPACEDIM==3      
        velxyz[0] = vel[0] * norm[0] + vel[1] * tan1[0] + vel[2] * tan2[0]; // x component
        velxyz[1] = vel[0] * norm[1] + vel[1] * tan1[1] + vel[2] * tan2[1]; // y component
        velxyz[2] = vel[0] * norm[2] + vel[1] * tan1[2] + vel[2] * tan2[2]; // z component
#else
        velxyz[0] = vel[0] * norm[0] + vel[1] * tan1[0];  // x component
        velxyz[1] = vel[0] * norm[1] + vel[1] * tan1[1];  // y component        
#endif

        Real P = prims[cls_t::QPRES];
        fluxw[cls_t::UMX] =  P*norm[0];
        fluxw[cls_t::UMY] =  P*norm[1];
#if AMREX_SPACEDIM==3          
        fluxw[cls_t::UMZ] =  P*norm[2];
#endif        
                                 
    }                                                           
                                                          
};
////////////////////////////////////////////////////////////////////////////
// adiabatic wall
template <typename cls_t>
class adiabatic_wall_t
{
  public:

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    adiabatic_wall_t() noexcept = default;	  
  // adiabatic_wall_t() {}
 
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    ~adiabatic_wall_t() noexcept = default;

    //~adiabatic_wall_t() {}

    static constexpr Real r43 = 4.0/3.0;
    static constexpr Real r13 = 1.0/3.0;

    // Eulerian flux
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void wall_flux (const auto &geomdata, int /*i*/, int /*j*/, int /*k*/, const Real norm[AMREX_SPACEDIM], 
      amrex::GpuArray<amrex::Real, cls_t::NPRIM>& prims, amrex::GpuArray<amrex::Real, cls_t::NCONS>& fluxw,const cls_t* cls) const noexcept {

      // printf(" oo Adiabatic wall \n ");
      // only set-non zero values (by default fluxw=0)                       
      Real P = prims[cls_t::QPRES];
      fluxw[cls_t::UMX] =  P*norm[0];
      fluxw[cls_t::UMY] =  P*norm[1];
#if AMREX_SPACEDIM==3          
      fluxw[cls_t::UMZ] =  P*norm[2];
#endif        
                                 
    }  
    // Viscous flux (stress and heat)  
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    void  wall_flux_diff(const auto &geomdata, int i, int j, int k, amrex::Real dis,const Real norm[AMREX_SPACEDIM],
      const Array4<Real>& q, amrex::GpuArray<amrex::Real, cls_t::NPRIM>& prims_w, 
      amrex::GpuArray<amrex::Real, cls_t::NCONS>& fluxw,const cls_t* cls, trans_parm_t const* ltransparm ) const noexcept  {
  
      // printf(" oo Adiabatic Viscous wall \n ");
      Real u[AMREX_SPACEDIM];
      u[0] = q(i,j,k,cls_t::QU); u[1] = q(i,j,k,cls_t::QV);
#if AMREX_SPACEDIM==3        
      u[2] = q(i,j,k,cls_t::QW);
#endif        
      // compute velocity derivatives normal direction
      Real dudn[AMREX_SPACEDIM]= {0.0};
      const Real dis_inv = 1.0/dis;
      for (int n = 0; n < AMREX_SPACEDIM; n++) {
        dudn[n] = u[n]*dis_inv;                
      }   

      // compute viscosity at the wall
#ifdef USE_PELEPHYSICS
      Real mu_w,cond_w,xi_w;
      const Real Tw = prims_w[cls_t::QT];
      // P,rho,Y
      //const Real P  = prims_w[cls_t::QPRES];
      const Real rho= prims_w[cls_t::QRHO];
      Real Y[NUM_SPECIES];
      for (int n = 0; n < NUM_SPECIES; ++n) { Y[n] = prims_w[cls_t::QFS + n]; }
      const bool get_xi = true, get_mu = true, get_lam = false;
      const bool get_Ddiag = false, get_chi = false;
      
      // pelePhysics       
      auto trans = pele::physics::PhysicsType::transport();   
      trans.transport(get_xi, get_mu, get_lam, get_Ddiag, get_chi, Tw, rho,
                     Y, nullptr, nullptr, mu_w, xi_w, cond_w, ltransparm);
      mu_w = mu_w*visc_cgs2si; // convert to SI units
#else
      Real mu_w  = cls->visc(prims_w[cls_t::QT]);
#endif
      // printf(" oo Adiabatic Viscous wall mu_w = %e  \n",mu_w);

      
      // coordinate transformation
#if AMREX_SPACEDIM==2      
      Real a1 = r43*norm[0]*norm[0] + norm[1]*norm[1];
      Real a2 = norm[0]*norm[0] + r43*norm[1]*norm[1]; 
      Real b1 = r13*norm[0]*norm[1];           
#elif AMREX_SPACEDIM==3        
      Real a1 = r43*norm[0]*norm[0]+ norm[1]*norm[1]  + norm[2]*norm[2];
      Real a2 = norm[0]*norm[0] + r43*norm[1]*norm[1] + norm[2]*norm[2];
      Real a3 = norm[0]*norm[0] + norm[1]*norm[1] + r43*norm[2]*norm[2];    
      Real b1 = r13*norm[0]*norm[1];             
      Real b2 = r13*norm[1]*norm[2]; 
      Real b3 = r13*norm[0]*norm[2];     
#endif        
       
      // adding stress and heat to flux 
#if AMREX_SPACEDIM==2
      fluxw[cls_t::UMX] -=  mu_w*(a1*dudn[0] + b1*dudn[1]);      
      fluxw[cls_t::UMY] -=  mu_w*(b1*dudn[0] + a2*dudn[1]);
#elif AMREX_SPACEDIM==3                  
      fluxw[cls_t::UMX] -=  mu_w*(a1*dudn[0] + b1*dudn[1] + b3*dudn[2]);
      fluxw[cls_t::UMY] -=  mu_w*(b1*dudn[0] + a2*dudn[1] + b2*dudn[2]);
      fluxw[cls_t::UMZ] -=  mu_w*(b3*dudn[0] + b2*dudn[1] + a3*dudn[2]);
#else
      amrex::Abort("Abort EBM only accepts 2D or 3D");  
      exit(1);              
#endif                
                                   
    }
                                                                  
};
////////////////////////////////////////////////////////////////////////////
// isothermal wall  ( Twall as input)
template < typename param, typename cls_t>
class isothermal_wall_t
{
  private:
  
    Real Twall = param::Twall; // wall temperature

  public:

    isothermal_wall_t() {}

    ~isothermal_wall_t() {}

    static constexpr Real r43 = 4.0/3.0;
    static constexpr Real r13 = 1.0/3.0;

    // Eulerian flux
    static void inline wall_flux(const auto &geomdata, int i, int j, int k, const Real norm[AMREX_SPACEDIM], 
      amrex::GpuArray<amrex::Real, cls_t::NPRIM>& prims, amrex::GpuArray<amrex::Real, cls_t::NCONS>& fluxw,const cls_t* cls) {      
      // only set-non zero values (by default fluxw=0)                       
      Real P = prims[cls_t::QPRES];
      fluxw[cls_t::UMX] =  P*norm[0];
      fluxw[cls_t::UMY] =  P*norm[1];
#if AMREX_SPACEDIM==3          
      fluxw[cls_t::UMZ] =  P*norm[2];
#endif        
 
      // species fluxes are 0

    }

    // Viscous flux (stress and heat)  
    static void inline wall_flux_diff(const auto &geomdata, int i, int j, int k, amrex::Real dis,const Real norm[AMREX_SPACEDIM],
      const Array4<Real>& q, amrex::GpuArray<amrex::Real, cls_t::NPRIM>& prims_w, 
      amrex::GpuArray<amrex::Real, cls_t::NCONS>& fluxw,const cls_t* cls) {

      // printf(" oo Isothermal Viscous wall \n ");              
      Real u[AMREX_SPACEDIM];
      u[0] = q(i,j,k,cls_t::QU); u[1] = q(i,j,k,cls_t::QV);
#if AMREX_SPACEDIM==3        
      u[2] = q(i,j,k,cls_t::QW);
#endif        
      // compute velocity derivatives normnal direction
      Real dudn[AMREX_SPACEDIM]= {0.0};
      const Real dis_inv = 1.0/dis;
      for (int n = 0; n < AMREX_SPACEDIM; n++) {
        dudn[n] = u[n]*dis_inv;                
      }   
      constexpr Real Tw = param::Twall;
     
      Real dTdn = (q(i,j,k,cls_t::QT) - Tw)*dis_inv;
      // compute viscosity and conductivity at the wall            
#ifdef USE_PELEPHYSICS
      Real mu_w,cond_w,xi_w;
      // P,rho,Y
      //const Real P  = prims_w[cls_t::QPRES];
      const Real rho= prims_w[cls_t::QRHO];
      Real Y[NUM_SPECIES];
      for (int n = 0; n < NUM_SPECIES; ++n) { Y[n] = prims_w[cls_t::QFS + n]; }
      const bool get_xi = true, get_mu = true, get_lam = true;
      const bool get_Ddiag = false, get_chi = false;

      // pelePhysics    
      auto trans = pele::physics::PhysicsType::transport();      
#if (PELEPVERSION==23)         
      trans_parms.allocate(); 
      auto const* ltransparm = trans_parms.device_trans_parm();
#else
      auto const* ltransparm = trans_parms.device_parm();
#endif    
    
      trans.transport(get_xi, get_mu, get_lam, get_Ddiag, get_chi, Tw, rho,
                     Y, nullptr, nullptr, mu_w, xi_w, cond_w, ltransparm);
      mu_w = mu_w*visc_cgs2si; // convert to SI units
      cond_w = cond_w*cond_cgs2si; // convert to SI units
     // xi_w = xi_w*visc_cgs2si; // convert to SI units                     
#else                     
      const Real mu_w   = cls->visc(Tw);
      const Real cond_w = cls->cond(Tw); 
#endif      
      // coordinate transformation
#if AMREX_SPACEDIM==2      
      Real a1 = r43*norm[0]*norm[0] + norm[1]*norm[1];
      Real a2 = norm[0]*norm[0] + r43*norm[1]*norm[1]; 
      Real b1 = r13*norm[0]*norm[1];           
#elif AMREX_SPACEDIM==3        
      Real a1 = r43*norm[0]*norm[0]+ norm[1]*norm[1]  + norm[2]*norm[2];
      Real a2 = norm[0]*norm[0] + r43*norm[1]*norm[1] + norm[2]*norm[2];
      Real a3 = norm[0]*norm[0] + norm[1]*norm[1] + r43*norm[2]*norm[2];     
      Real b1 = r13*norm[0]*norm[1];            
      Real b2 = r13*norm[1]*norm[2]; 
      Real b3 = r13*norm[0]*norm[2];     
#endif        
       
      // adding stress and heat to flux 
#if AMREX_SPACEDIM==2
      fluxw[cls_t::UMX] -=  mu_w*(a1*dudn[0] + b1*dudn[1]);      
      fluxw[cls_t::UMY] -=  mu_w*(b1*dudn[0] + a2*dudn[1]);
#elif AMREX_SPACEDIM==3                  
      fluxw[cls_t::UMX] -=  mu_w*(a1*dudn[0] + b1*dudn[1] + b3*dudn[2]);
      fluxw[cls_t::UMY] -=  mu_w*(b1*dudn[0] + a2*dudn[1] + b2*dudn[2]);
      fluxw[cls_t::UMZ] -=  mu_w*(b3*dudn[0] + b2*dudn[1] + a3*dudn[2]);
#else
      amrex::Abort("Abort EBM only accepts 2D or 3D");  
      exit(1);              
#endif                
      fluxw[cls_t::UET] -=  cond_w*dTdn;

      // all species fluxes are 0 ! 
        
    }    

};
////////////////////////////////////////////////////////////////////////////
// user wall

#endif

