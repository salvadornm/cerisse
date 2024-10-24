#ifndef Transport2_H_
#define Transport2_H_

////////////////////////////////TRANSPORT/////////////////////////////////
template <typename param>
class transport_const_t {
  private:
  Real visc_ref = param::viscosity;
  Real cond_ref = param::conductivity;
  Real xi_ref   = 0.0;

  public:

  AMREX_GPU_DEVICE AMREX_FORCE_INLINE Real visc(const Real& T) const {
    return visc_ref;
  } 
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE Real cond(Real& T) const {
    return cond_ref;
  } 
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE Real xi(Real& T) const {
    return xi_ref;
  } 
};

class transport_suth_t {

  private:
  // Sutherland's fit from Computational Fluid Mechanics and Heat Transfer
  Real visc_ref   = 1.458e-6; //viscosity parameters
  Real Tvisc_ref  = 110.4;
  Real cond_ref   = 2.495e-3; //conductivity parameters
  Real Tcond_ref  = 194.0;
  Real xi_ref     = 0.0; // bulk

  public:
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE Real visc(const Real& T) const {
    return visc_ref * T * sqrt(T) / (Tvisc_ref + T);
  }

  AMREX_GPU_DEVICE AMREX_FORCE_INLINE Real cond(Real& T) const {
    return cond_ref * T * sqrt(T) / (Tcond_ref + T);
  }

  AMREX_GPU_DEVICE AMREX_FORCE_INLINE Real xi(Real& T) const {
    return xi_ref;
  }

};

////////////////////////////////////////////////////////////////////////////////

#endif