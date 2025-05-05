//--------------------------------------------------------------------------//
// param is a struct withteh follwoing options

// isothermal wall  ( Twall as input)
template < typename param, typename cls_t>
class ibm_isothermal_slip_wall_t
{
  private:
  
    Real Twall = param::Twall; // wall temperature

  public:

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void computeIB(Array2D<Real,0,eorder_tparm+1,0,cls_t::NPRIM-1>& q, const cls_t* cls) {
    
    // slip velocity (in local coordinates)
    q(1,cls_t::QU) = 0.0_rt; // un
    q(1,cls_t::QV) = q(2,cls_t::QV); // ut1
    q(1,cls_t::QW)  = q(2,cls_t::QW); // ut2

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
// adiabatic wall  ( Twall as input)
template < typename param, typename cls_t>
class ibm_adiabatic_slip_wall_t
{
  private:
  
  public:

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void computeIB(Array2D<Real,0,eorder_tparm+1,0,cls_t::NPRIM-1>& q, const cls_t* cls) {
    
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

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void computeIB(Array2D<Real,0,eorder_tparm+1,0,cls_t::NPRIM-1>& q, const cls_t* cls) {
        
    for (int n = 0; n <= cls_t::QLS; ++n) {
      q(1,n) = param::alpha[n]*q(2,n) + param::beta[n]; 
    }      
    
  }
};    




