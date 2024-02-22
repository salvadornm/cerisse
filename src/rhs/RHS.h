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

void inline static compute() {
  std::cout << "compute_rhs" << std::endl;
  exit(0);
};

// void CNS::compute_rhs(MultiFab& statemf, Real dt, FluxRegister* fr_as_crse, FluxRegister* fr_as_fine) {
//   BL_PROFILE("CNS::compute_rhs()");

//   // Variables
//   // TODO: introduce a struct for these variables?
//   const PROB::ProbClosures& cls_d = *CNS::d_prob_closures;
//   const PROB::ProbClosures& cls_h = *CNS::h_prob_closures;
//   const PROB::ProbParm& parms = *d_prob_parm;

//   for (MFIter mfi(statemf, false); mfi.isValid(); ++mfi) {
//     Array4<Real> const& state = statemf.array(mfi);

//     const Box& bxgnodal = mfi.grownnodaltilebox(-1, 0);  // extent is 0,N_cell+1
//     const Box& bxg = mfi.growntilebox(cls_d.NGHOST);

//     FArrayBox primf(bxg, cls_d.NPRIM , The_Async_Arena());
//     FArrayBox tempf(bxgnodal, cls_d.NCONS, The_Async_Arena());
//     Array4<Real> const& temp = tempf.array();
//     Array4<Real> const& prims= primf.array();

//     // We want to minimise function calls. So, we call prims2cons, flux and
//     // source term evaluations once per fab from CPU, to be run on GPU.
//     cls_h.cons2prims(mfi, state, prims);

//     // Fluxes including boundary/discontinuity corrections
//     // Note: we are over-writing state (cons) with flux derivative
//     prob_rhs.eflux(geom, mfi, prims, temp, state, cls_h);
//     prob_rhs.dflux(); //(prims,cons,nflx)

//     // Source terms, including update mask (e.g inside IB)
//     // prob_rhs.src(prims,cons,nflx,rhs)
//   }
// }

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