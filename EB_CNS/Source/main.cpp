#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Amr.H>

#if (AMREX_SPACEDIM > 1) //1D cannot have EB
#include <AMReX_EB2.H>
#endif

#include "CNS.H"

using namespace amrex;

amrex::LevelBld* getLevelBld ();

#if (AMREX_SPACEDIM > 1) //1D cannot have EB
void initialize_EB2 (const Geometry& geom, const int required_level, const int max_level);
#endif

int main (int argc, char* argv[])
{
  amrex::Initialize(argc,argv);

  BL_PROFILE_VAR("main()", pmain);

  double timer_tot = amrex::second();
  double timer_init = 0.;
  double timer_advance = 0.;

  int  max_step;
  Real strt_time;
  Real stop_time;

  {
    ParmParse pp;

    max_step  = -1;
    strt_time =  0.0;
    stop_time = -1.0;

    pp.query("max_step",max_step);
    pp.query("strt_time",strt_time);
    pp.query("stop_time",stop_time);
  }

  if (strt_time < 0.0) {
    amrex::Abort("MUST SPECIFY a non-negative strt_time");
  }

  if (max_step < 0 && stop_time < 0.0) {
    amrex::Abort("Exiting because neither max_step nor stop_time is non-negative.");
  }

  {
    timer_init = amrex::second();

    Amr amr(getLevelBld());

#if (AMREX_SPACEDIM > 1) //1D cannot have EB
    AmrLevel::SetEBSupportLevel(EBSupport::full);// need both area and volume fractions
    AmrLevel::SetEBMaxGrowCells(NUM_GROW, 4, 2); // 5 for ebcellflags, 4 for vfrac, 2 is not used for EBSupport::volume
    initialize_EB2(amr.Geom(amr.maxLevel()), amr.maxLevel(), amr.maxLevel());
#endif

    amr.init(strt_time,stop_time);

    timer_init = amrex::second() - timer_init;

    timer_advance = amrex::second();

    while ((amr.okToContinue()) && 
           (amr.levelSteps(0) < max_step || max_step < 0) &&
           (amr.cumTime() < stop_time || stop_time < 0.0))
    {
      // Do a coarse timestep. Recursively calls timeStep()
      amr.coarseTimeStep(stop_time);
    }

    timer_advance = amrex::second() - timer_advance;

    // Write final checkpoint and plotfile
    if (amr.stepOfLastCheckPoint() < amr.levelSteps(0)) {
      amr.checkPoint();
    }

    if (amr.stepOfLastPlotFile() < amr.levelSteps(0)) {
      amr.writePlotFile();
    }
  }

  timer_tot = amrex::second() - timer_tot;

  ParallelDescriptor::ReduceRealMax({timer_tot, timer_init, timer_advance},
                                    ParallelDescriptor::IOProcessorNumber());

  amrex::Print() << "Run Time total   = " << timer_tot     << "\n"
                 << "Run Time init    = " << timer_init    << "\n"
                 << "Run Time advance = " << timer_advance << "\n";

  BL_PROFILE_VAR_STOP(pmain);

  amrex::Finalize();

  return 0;
}