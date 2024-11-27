#ifndef RHS_H_
#define RHS_H_

//#include <Index.h>

// Euler numerical methods
#include <Weno.h>
#include <CentralKEEP.h>
#include <Riemann.h>
#include <Rusanov.h>
#include <Skew.h>

// viscous templates
#include <DiffusionCD.h>
#include <viscous.h>

#ifdef USE_PELEPHYSICS
#include "react.h"
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
  void eflux_ibm(Args&&... args){}
#else  
  void eflux(Args&&... args){}
#endif  
};

// no diffusive flux
class no_diffusive_t
{
public:
  template<typename... Args>
#if (AMREX_USE_GPIBM || CNS_USE_EB )   
  void dflux_ibm(Args&&... args) {}
#else
  void dflux(Args&&... args) {}  
#endif  
};

// no source
class no_source_t
{
public:
  template<typename... Args>
  void src(Args&&... args) {}
};

#endif