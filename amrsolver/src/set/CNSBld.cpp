#include <AMReX_LevelBld.H>
#include <CNS.h>

using namespace amrex;

// LevelBls is a virtual class.
class CNSBld : public LevelBld {
  void variableSetUp() override;
  void variableCleanUp() override;
  AmrLevel* operator()() override;
  AmrLevel* operator()(Amr& papa, int lev, const Geometry& level_geom,
                       const BoxArray& ba, const DistributionMapping& dm,
                       Real time) override;
};

CNSBld CNS_bld;

LevelBld* getLevelBld() { return &CNS_bld; }

void CNSBld::variableSetUp() { CNS::variableSetUp(); }

void CNSBld::variableCleanUp() { CNS::variableCleanUp(); }

AmrLevel* CNSBld::operator()() { return new CNS; }

// This function overwrites a virtual function in LevelBld
AmrLevel* CNSBld::operator()(Amr& papa, int lev, const Geometry& level_geom,
                             const BoxArray& ba, const DistributionMapping& dm,
                             Real time) {
  return new CNS(papa, lev, level_geom, ba, dm, time);
}
