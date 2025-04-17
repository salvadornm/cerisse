#ifndef walltypes_H_
#define walltypes_H_

#include <AMReX_EB2.H>
#include <AMReX_EB2_IF.H>
#include <AMReX_ParmParse.H>

#include <algorithm>
#include <cmath>

#include <EBMultiFab.h>


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

    adiabatic_wall_t() {}

    ~adiabatic_wall_t() {}

    // Eulerian flux
    static void inline wall_flux(const auto &geomdata, int i, int j, int k, const Real norm[AMREX_SPACEDIM], 
      amrex::GpuArray<amrex::Real, cls_t::NPRIM>& prims, amrex::GpuArray<amrex::Real, cls_t::NCONS>& fluxw,const cls_t* cls) {

        // printf(" oo Adiabatic wall \n ");
        // only set-non zero values (by default fluxw=0)                       
        Real P = prims[cls_t::QPRES];
        fluxw[cls_t::UMX] =  P*norm[0];
        fluxw[cls_t::UMY] =  P*norm[1];
#if AMREX_SPACEDIM==3          
        fluxw[cls_t::UMZ] =  P*norm[2];
#endif        
                                 
    }  
    // Viscous flux (stress and heat)  (tau) 
    static void inline wall_flux_diff(const auto &geomdata, int i, int j, int k, const Real norm[AMREX_SPACEDIM],
        const Array4<Real>& q, amrex::GpuArray<amrex::Real, cls_t::NPRIM>& prims_w, 
        amrex::GpuArray<amrex::Real, cls_t::NCONS>& fluxw,const cls_t* cls) {
  
        // printf(" oo Adiabatic Viscous wall \n ");
        // tangential vectors 
        Real tan1[AMREX_SPACEDIM];
        tan1[0]  =  norm[1];tan1[1]  = -norm[0];
#if AMREX_SPACEDIM==3        
        tan1[2] = 0.0;
        Real tan2[AMREX_SPACEDIM];
        tan2[0] = tan1[1] * norm[2] - tan1[2] * norm[1]; // x component
        tan2[1] = tan1[2] * norm[0] - tan1[0] * norm[2]; // y component
        tan2[2] = tan1[0] * norm[1] - tan1[1] * norm[0]; // z component      
#endif        
        Real u[AMREX_SPACEDIM];
        u[0] = q(i,j,k,cls_t::QU); u[1] = q(i,j,k,cls_t::QV);
#if AMREX_SPACEDIM==3        
        u[2] = q(i,j,k,cls_t::QW);
#endif        
        Real dudn,dvdn,dwdn,dTdn;

        // compute conductivity and viscosity
        Real mu_  = cls->visc(prims_w[cls_t::QT]);
        Real lam  = cls->cond(prims_w[cls_t::QT]);  

        fluxw[cls_t::UMX] +=  0.0;
        fluxw[cls_t::UMY] +=  0.0;
  #if AMREX_SPACEDIM==3          
        fluxw[cls_t::UMZ] +=  0.0;
  #endif        
        fluxw[cls_t::UET] +=  lam*dTdn;

                                   
      }
    
    
                                                          
};
////////////////////////////////////////////////////////////////////////////
// isothermal wall  (requires Twall)

////////////////////////////////////////////////////////////////////////////
// user wall

#endif

