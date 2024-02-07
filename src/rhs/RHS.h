#ifndef RHS_H_
#define RHS_H_

#include <Index.h>
#include <Weno.h>
#include <CentralKEEP.h>
#include <Riemann.h>

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
  void eflux(){};
};

// no diffusive flux
class no_diffusive_t
{
public:
  void dflux() {}
};

// no source
class no_source_t
{
public:
  void src() {}
};

#endif