#ifndef TransportPele_H_
#define TransportPele_H_

#ifdef USE_PELEPHYSICS
#include <PelePhysics.H>

static pele::physics::transport::TransportParams<
    pele::physics::PhysicsType::transport_type> trans_parms;

#endif
////////////////////////////////TRANSPORT/////////////////////////////////

class transport_Pele_t {

  private:
  
  public:

  // constructor
  AMREX_GPU_HOST_DEVICE
  transport_Pele_t()
  {
   // trans_parms.allocate(); 
  }


  // These are dummy calls, Cerisse not expected to call these functions directly using PelePhsyics
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE Real visc(const Real& T) const {
    return -99;
  }

  AMREX_GPU_DEVICE AMREX_FORCE_INLINE Real cond(Real& T) const {
    return -99;
  }

  AMREX_GPU_DEVICE AMREX_FORCE_INLINE Real xi(Real& T) const {
    return -99;
  }

#ifdef USE_PELEPHYSICS

  //static pele::physics::transport::TransportParams<
  //  pele::physics::PhysicsType::transport_type> trans_parms;
    
#endif

};

////////////////////////////////////////////////////////////////////////////////

#endif