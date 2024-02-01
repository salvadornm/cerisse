#ifndef Transport_H_
#define Transport_H_

////////////////////////////////TRANSPORT/////////////////////////////////
////// Viscosity
class visc_const_t {
 private:
 public:
  Real visc_ref = 1.458e-6;  // Viscosity reference value

  AMREX_GPU_DEVICE AMREX_FORCE_INLINE Real visc(const Real& T) const {
    return visc_ref;
  }
};

class visc_suth_t {
 private:
  // Sutherland's fit from Computational Fluid Mechanics and Heat Transfer
  Real visc_ref = 1.458e-6;
  Real Tvisc_ref = 110.4;

 public:
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE Real visc(const Real& T) const {
    return visc_ref * T * sqrt(T) / (Tvisc_ref + T);
  }
};

////// Conductivity
class cond_const_t {
 private:
 public:
  Real cond_ref = 1.458e-6;

  AMREX_GPU_DEVICE AMREX_FORCE_INLINE Real cond(Real& T) const {
    return cond_ref;
  }

#ifdef AMREX_USE_GPU
  AMREX_FORCE_INLINE Real cond_cpu(Real& T) const { return cond_ref; }
#endif
};

// class cond_const_pr_t private: ThermodynamicsBase {
//   private:
//   public:
//     Real visc_ref = 1.458e-6; // Viscosity reference value

//     AMREX_GPU_DEVICE AMREX_FORCE_INLINE Real cond(Real& T) const {
//       // cp*;
//     return visc_ref;}
// };

class cond_suth_t {
 private:
  // Sutherland's fit from Computational Fluid Mechanics and Heat Transfer
  Real cond_ref = 2.495e-3;
  Real Tcond_ref = 194.0;

 public:
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE Real cond(Real& T) const {
    return cond_ref * T * sqrt(T) / (Tcond_ref + T);
  }

#ifdef AMREX_USE_GPU
  AMREX_FORCE_INLINE Real cond_cpu(Real& T) const {
    return cond_ref * T * sqrt(T) / (Tcond_ref + T);
  }
#endif
};
////////////////////////////////////////////////////////////////////////////////

#endif