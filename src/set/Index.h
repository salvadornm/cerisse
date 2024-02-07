#ifndef INDEX_H_
#define INDEX_H_

#include <AMReX_GpuContainers.H>

struct indicies_t {

public:
// if statement based selection on types of indicies
// Independent (solved) variables

// Static constexpr not needed, but this allows compile time use of indicies if needed.
static constexpr int URHO=0;
static constexpr int UMX=1;
static constexpr int UMY=2;
static constexpr int UMZ=3;
static constexpr int UET=4;

// Dependent (derived) variables
static constexpr int QRHO=0;
static constexpr int QU=1;
static constexpr int QV=2;
static constexpr int QW=3;
static constexpr int QT=4;
static constexpr int QPRES=5;

// Vector<std::string> prim_vars_names={"Pressure","Temperature","Xvel","Yvel","Zvel"};
// amrex::Gpu::ManagedVector<std::string> cons_vars_names={"Density","Xmom","Ymom","Zmom","Energy"};

// amrex::Gpu::ManagedVector<char[2]> cons_vars_names;
// ={"Density","Xmom","Ymom","Zmom","Energy"};
      // pr        ux      uy     uz      et
// amrex::Gpu::ManagedVector<int> cons_vars_type;
//  ={0,1,2,3,0};
// amrex::ManagedVector<std::string> prim_vars_names={"Pressure","Temperature"};

static constexpr int NCONS=5;
static constexpr int NPRIM=6;
static constexpr int NGHOST=3; // TODO: make it an automatic parameter
};

#endif