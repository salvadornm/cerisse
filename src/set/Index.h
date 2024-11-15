#ifndef INDEX_H_
#define INDEX_H_

#include <AMReX_GpuContainers.H>

#ifdef USE_PELEPHYSICS
#include <PelePhysics.H>
#include <mechanism.H>
#else
static constexpr int NUM_SPECIES = 1;
#endif

struct indicies_t {

public:
// if statement based selection on types of indicies
// Independent (solved) variables

// Static constexpr not needed, but this allows compile time use of indicies if needed.
static constexpr int UMX=0;
static constexpr int UMY=1;
static constexpr int UMZ=2;
static constexpr int UET=3;
static constexpr int URHO=4;
static constexpr int UFS=URHO;

static amrex::Vector<std::string> get_cons_vars_names() {  
#ifdef USE_PELEPHYSICS
  amrex::Vector<std::string> spec_names;
  pele::physics::eos::speciesNames<pele::physics::PhysicsType::eos_type>(spec_names);
  amrex::Vector<std::string> names = {"Xmom", "Ymom", "Zmom", "Energy"};  
  for (int i = 0; i < NUM_SPECIES; ++i) {
    names.push_back("rho_" + spec_names[i]);
  }
  return names;
#else
  return {"Xmom", "Ymom", "Zmom", "Energy", "Density"};
#endif
}

static amrex::Vector<int> get_cons_vars_type() {
  amrex::Vector<int> types = {1, 2, 3, 0};
  for (int i = 0; i < NUM_SPECIES; ++i) {
    types.push_back(0);
  }
  return types;
}

// Dependent (derived) variables
static constexpr int QRHO=0;
static constexpr int QU=1;
static constexpr int QV=2;
static constexpr int QW=3;
static constexpr int QT=4;
static constexpr int QPRES=5;
static constexpr int QC=6;
static constexpr int QG=7;
static constexpr int QEINT=8;
static constexpr int QFS=9;

// Transport coefficients
static constexpr int CMU=0;   // dynamic viscosity
static constexpr int CXI=1;   // bulk viscosity
static constexpr int CLAM=2;  // heat conductivity 
static constexpr int CRHOD=3; // species diffusivity (rho*D)

static constexpr int NCONS=UFS + NUM_SPECIES;
static constexpr int NPRIM=QFS + NUM_SPECIES;
static constexpr int NWAVES=2 + NUM_SPECIES;
static constexpr int NCOEF= 3 + NUM_SPECIES;
static constexpr int NGHOST=3; // TODO: make it an automatic parameter

// Statistic space empty  
static constexpr int NSTAT = 0;    
static constexpr int NSTAT_VEL = 0;    
static constexpr int NSTAT_THERM = 0;    



};

#endif