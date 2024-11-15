#ifndef CLOSURES_H_
#define CLOSURES_H_

#include <Index.h>
#include <Index_stat.h>

#include <NumParam.h>
#include <Thermodynamics.h>
#include <Transport.h>
#include <Transport2.h>

//#include <Turbulence.h>

using namespace amrex;


// template <typename idx, typename Thermo, typename MolTransport, typename... others>
// class closures_dt : public idx, public Thermo,  public MolTransport, public others... {

template <typename idx, typename Thermo, typename... others>
class closures_dt : public idx, public Thermo, public others... {

//template <typename idx, typename Visc, typename Cond, typename Thermo, typename... others>
//class closures_dt : public idx, public Cond, public Visc, public Thermo,  public others... {
 private:
 public:

};
#endif