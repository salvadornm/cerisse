#ifndef RHS_H_
#define RHS_H_

//#include <Index.h>
#include <CNS.h>

// Euler numerical methods
#include <Weno.h>
#include <CentralKEEP.h>
#include <CentralDif.h>
#include <Riemann.h>
#include <Rusanov.h>
#include <Skew.h>

// viscous templates
#include <DiffusionCD.h>
#include <viscous.h>
#include <viscousLES.h>


#ifdef USE_PELEPHYSICS
#include "react.h"
#include "react_source.h"
#include "react_sourceLES.h"
#endif

// _dt stands for derived type
template <typename euler, typename diffusive, typename source>
class rhs_dt : public euler, public diffusive, public source
{
private:
public:
};

// no euler flux
class no_euler_t
{
public:
  template<typename... Args>
#if (AMREX_USE_GPIBM || CNS_USE_EB )  
  // void eflux_ibm(Args&&... args){}
  void eflux_ibm(const Geometry& /*geom*/, const MFIter& /*mfi*/,
                    const Array4<Real>& /*prims*/, std::array<FArrayBox*, AMREX_SPACEDIM> const /*&flxt*/,
                    const Array4<Real>& /*cons*/, Args&&... args ) { }
#else  
  //void eflux(Args&&... args){}
  void eflux(const Geometry& /*geom*/, const MFIter& /*mfi*/,
            const Array4<Real>& /*prims*/, std::array<FArrayBox*, AMREX_SPACEDIM> const /*&flxt*/,            
            const Array4<Real>& /*rhs*/, Args&&... args) { }
#endif  
};

// no diffusive flux
class no_diffusive_t
{
public:

  // No-op init, so ProbRHS::init_coeffs() is always valid
  AMREX_GPU_HOST
  void init_coeffs() {}

  template<typename... Args>
#if (AMREX_USE_GPIBM || CNS_USE_EB )   
  //void dflux_ibm(Args&&... args) {}
  void dflux_ibm(const Geometry& geom, const MFIter& mfi,
            const Array4<Real>& prims, std::array<FArrayBox*, AMREX_SPACEDIM> const &flxt,            
            const Array4<Real>& rhs, Args&&... args) { }

#else
  //void dflux(Args&&... args) {}  
  void dflux(const Geometry& geom, const MFIter& mfi,
            const Array4<Real>& prims, std::array<FArrayBox*, AMREX_SPACEDIM> const &flxt,            
            const Array4<Real>& rhs, Args&&... args) { }

#endif  

  // no-op RZ geometric viscous source (no diffusion => no hoop stress)
  void inline rz_geometric_source(const Geometry& /*geom*/, const MFIter& /*mfi*/,
            const Array4<Real>& /*prims*/, const Array4<Real>& /*state*/,
            const auto* /*cls*/) { }
};

// no source
class no_source_t
{
public:
  template<typename... Args>
  void src(Args&&... args) {}
};

#endif