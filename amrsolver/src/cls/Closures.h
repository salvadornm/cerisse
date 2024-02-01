#ifndef CLOSURES_H_
#define CLOSURES_H_

#include <Index.h>
#include <Thermodynamics.h>
#include <Transport.h>

using namespace amrex;

template <typename idx, typename Visc, typename Cond, typename Thermo, typename... others>
class closures_dt : public idx, public Cond, public Visc, public Thermo,  public others... {
 private:
 public:

};
#endif