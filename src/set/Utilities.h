#ifndef UTILITY_H_
#define UTILITY_H_

#include <AMReX_AmrLevel.H>
#include <AMReX_FluxRegister.H>
#include <prob.h>
#include <CNSconstants.h>

#ifdef USE_PELEPHYSICS
#include <PMFData.H>
#include <PMF.H>
#endif

using namespace amrex;

// interface between PelePhysics<-> Cerisse as well as other utilities

class Utility{

  public:
  
  Utility()
  {};

  ~Utility()
  {};

  void test(){
    std::cout << " ***  Utility test ***" << std::endl;
  }

// PMF: Premixed Flame Initilialization
#ifdef USE_PELEPHYSICS

  pele::physics::PMF::PmfData pmfData;

  // initialize PMF data
  void initPMF() {
    pmfData.initialize();
  }
#endif


  // Turbulence



  private:

};

#endif  
 

