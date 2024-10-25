#ifndef RHS_H_
#define RHS_H_

#include <Index.h>

// numerical methods
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
// later on change this from rhs_dt to cns_dt.
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
  void eflux(Args&&... args){}
  template<typename... Args>
  void eflux_ibm(Args&&... args){}
};

// no diffusive flux
class no_diffusive_t
{
public:
  template<typename... Args>
  void dflux(Args&&... args) {}
  template<typename... Args>
  void dflux_ibm(Args&&... args) {}
};

// no source
class no_source_t
{
public:
  template<typename... Args>
  void src(Args&&... args) {}
};

#endif