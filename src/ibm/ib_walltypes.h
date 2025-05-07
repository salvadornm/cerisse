#ifndef ib_walltypes_H_
#define ib_walltypes_H_

#include <IBMultiFab.h>
#include <AMReX_GpuContainers.H>
#include <AMReX_IntVect.H>
#include <AMReX_StateDescriptor.H>
#include <AMReX_Derive.H>

//--------------------------------------------------------------------------//
// \brief templates for different wall types for IB 
// param is a struct with the follwoing options
// \param Twall : wall temperture (required for isothermal)
// \param alpha : array of coefficients (reuired for generic bc)
// \param beta  : array of coefficients (reuired for generic bc)
// the compute IB calculates surface primitive variables:
//        normal and tangential velocities, P, T and mass fractions
//        as a function of x,y,z, normal and interpolated vars
//--------------------------------------------------------------------------//

// isothermal slip wall
template < typename param, typename cls_t>
class ibm_isothermal_slip_wall_t
{
  private:

  public:
  
  static constexpr Real Twall   = param::Twall;
  static const int eorder_tparm = param::extrap_order;

  
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  static void compute_surfIB(const Array1D<Real,0,AMREX_SPACEDIM-1>& /*xyz*/,const Array1D<Real,0,AMREX_SPACEDIM-1>& /*norm*/,
                      Array2D<Real,0,eorder_tparm+1,0,cls_t::NPRIM-1>& q, const cls_t* cls) {
    
    // slip velocity (in local coordinates)
    q(1,cls_t::QU) = 0.0_rt; // un
    q(1,cls_t::QV) = q(2,cls_t::QV); // ut1
    q(1,cls_t::QW) = q(2,cls_t::QW); // ut2

    Real Yw[NUM_SPECIES]={0.0};

    // zerograd pressure    
    q(1,cls_t::QPRES) = q(2,cls_t::QPRES); 
    // wall temperature
    q(1,cls_t::QT)    = param::Twall;

#if NUM_SPECIES > 1    
    Real sumY = 0.0;
    for (int n = 0; n < NUM_SPECIES; ++n) {
      Yw[n]   =  q(2,cls_t::QFS+n);
      sumY + = sumY;
    }
    for (int n = 0; n < NUM_SPECIES; ++n) { 
      q(1,cls_t::QFS+n)   =  Yw[n]/sumY;
    }
#endif                      
  }
};    
//--------------------------------------------------------------------------//
// isothermal non-slip wall
template < typename param, typename cls_t>
class ibm_isothermal_nonslip_wall_t
{
  private:
  
  public:

    static constexpr Real Twall = param::Twall;
    static const int eorder_tparm = param::extrap_order;

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static void compute_surfIB(const Array1D<Real,0,AMREX_SPACEDIM-1>& /*xyz*/,const Array1D<Real,0,AMREX_SPACEDIM-1>& /*norm*/,
      Array2D<Real,0,eorder_tparm+1,0,cls_t::NPRIM-1>& q, const cls_t* cls) {

      // slip velocity (in local coordinates)
      q(1,cls_t::QU) = 0.0; // un
      q(1,cls_t::QV) = 0.0; // ut1
      q(1,cls_t::QW) = 0.0; // ut2

      Real Yw[NUM_SPECIES]={0.0};

      // zerograd pressure    
      q(1,cls_t::QPRES) = q(2,cls_t::QPRES); 
      // wall temperature
      q(1,cls_t::QT)    = param::Twall;

#if NUM_SPECIES > 1    
      Real sumY = 0.0;
      for (int n = 0; n < NUM_SPECIES; ++n) {
        Yw[n]   =  q(2,cls_t::QFS+n);
        sumY + = sumY;
      }
      for (int n = 0; n < NUM_SPECIES; ++n) { 
        q(1,cls_t::QFS+n)   =  Yw[n]/sumY;
      }
#endif                      
    }
};    
//--------------------------------------------------------------------------//
// adiabatic slip wall  
template < typename param, typename cls_t>
class ibm_adiabatic_slip_wall_t
{
  private:
  
  public:

    static const int eorder_tparm = param::extrap_order;

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static void compute_surfIB(const Array1D<Real,0,AMREX_SPACEDIM-1>& /*xyz*/,const Array1D<Real,0,AMREX_SPACEDIM-1>& /*norm*/,
      Array2D<Real,0,eorder_tparm+1,0,cls_t::NPRIM-1>& q, const cls_t* cls) {

      // slip velocity (in local coordinates)
      q(1,cls_t::QU) = 0.0_rt; // un
      q(1,cls_t::QV) = q(2,cls_t::QV); // ut1
      q(1,cls_t::QW)  = q(2,cls_t::QW); // ut2

      Real Yw[NUM_SPECIES]={0.0};

      // zerograd pressure    
      q(1,cls_t::QPRES) = q(2,cls_t::QPRES); 
      // wall temperature
      q(1,cls_t::QT)    = q(2,cls_t::QT);

#if NUM_SPECIES > 1    
      Real sumY = 0.0;
      for (int n = 0; n < NUM_SPECIES; ++n) {
        Yw[n]   =  q(2,cls_t::QFS+n);
        sumY + = sumY;
      }
      for (int n = 0; n < NUM_SPECIES; ++n) { 
        q(1,cls_t::QFS+n)   =  Yw[n]/sumY;
      }
#endif                      
    }
};    
//--------------------------------------------------------------------------//
// general boundary
// imposes bc of teh form phi(1) =alpha*phi(2) + beta
// alpha = 1  beta=0      dphi/dn=0   Neumann
// alpha = 0  beta=PHIBC  phi =PHIBC  Dirichlet

template < typename param, typename cls_t>
class ibm_general_wall_t
{
  private:

  public:

    static const int eorder_tparm = param::extrap_order;

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static void compute_surfIB(const Array1D<Real,0,AMREX_SPACEDIM-1>& /*xyz*/,const Array1D<Real,0,AMREX_SPACEDIM-1>& /*norm*/,
      Array2D<Real,0,eorder_tparm+1,0,cls_t::NPRIM-1>& q, const cls_t* cls) {

      for (int n = 0; n <= cls_t::QLS; ++n) {
       q(1,n) = param::alpha[n]*q(2,n) + param::beta[n]; 
      }      
    }
};    

//-----------------------------------------------------------------------------//
#endif

