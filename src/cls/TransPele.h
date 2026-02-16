#ifndef TRANSPORTPELE_H_
#define TRANSPORTPELE_H_

#ifdef USE_PELEPHYSICS
#include <PelePhysics.H>


#if (PELEPVERSION==23)
// v23
static  pele::physics::transport::TransportParams<
        pele::physics::PhysicsType::transport_type> trans_parms;

// check (not validated v23)
using trans_parm_t =
 	pele::physics::transport::TransParm<
	pele::physics::PhysicsType::eos_type,
	pele::physics::PhysicsType::transport_type>;	

#else
// v25       
extern  pele::physics::PeleParams<
        pele::physics::transport::TransParm<
        pele::physics::PhysicsType::eos_type,
        pele::physics::PhysicsType::transport_type > > trans_parms;

using trans_parm_t =
	typename pele::physics::transport::TransParm<
	pele::physics::PhysicsType::eos_type,
	pele::physics::PhysicsType::transport_type>;

#endif



#endif
////////////////////////////////TRANSPORT/////////////////////////////////

class transport_Pele_t {

  private:

  // default values (not to be used)
  Real visc_ref = 1.458e-6;
  Real cond_ref = 2.495e-3;
  Real xi_ref = 0.0;
  
  public:

  // constructor
  AMREX_GPU_HOST_DEVICE
  transport_Pele_t()
  {
   // trans_parms.allocate();    
    //trans_parms.initialize();
  }


  // These are dummy calls, Cerisse not expected to call these functions directly using PelePhsyics
  // WARNING !! This will cause a problem if you use EBM+PelePhysics
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE Real visc(const Real& T) const {
    return (visc_ref);
  }

  AMREX_GPU_DEVICE AMREX_FORCE_INLINE Real cond(const Real& T) const {
    return (cond_ref);
  }

  AMREX_GPU_DEVICE AMREX_FORCE_INLINE Real xi(const Real& T) const {
    return (xi_ref);
  }

#ifdef USE_PELEPHYSICS

  //static pele::physics::transport::TransportParams<
  //  pele::physics::PhysicsType::transport_type> trans_parms;
    
#endif

};

////////////////////////////////////////////////////////////////////////////////

#endif  // TRANSPORTPELE_H_
