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

// Dependent (derived) variables
static constexpr int QRHO=0;
static constexpr int QU=1;
static constexpr int QV=2;
static constexpr int QW=3;
static constexpr int QT=4;
static constexpr int QPRES=5;
static constexpr int QC=6;
static constexpr int QG=7;
static constexpr int QFS=8;

// Vector<std::string> prim_vars_names={"Pressure","Temperature","Xvel","Yvel","Zvel"};
// amrex::Gpu::ManagedVector<std::string> cons_vars_names={"Density","Xmom","Ymom","Zmom","Energy"};

// amrex::Gpu::ManagedVector<char[2]> cons_vars_names;
// ={"Density","Xmom","Ymom","Zmom","Energy"};
      // pr        ux      uy     uz      et
// amrex::Gpu::ManagedVector<int> cons_vars_type;
//  ={0,1,2,3,0};
// amrex::ManagedVector<std::string> prim_vars_names={"Pressure","Temperature"};

static constexpr int NCONS=UFS + NUM_SPECIES;
static constexpr int NPRIM=QFS + NUM_SPECIES;
static constexpr int NGHOST=2; // TODO: make it an automatic parameter
};

#endif