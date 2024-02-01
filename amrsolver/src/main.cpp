#include <AMReX.H>
#include <AMReX_Amr.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <CNS.h>
#include <prob.h>

#ifdef AMREX_USE_GPIBM
#include <IBM.h>
#endif

using namespace amrex;

amrex::LevelBld* getLevelBld();

int main(int argc, char* argv[]) {
  BL_PROFILE_VAR("main()", pmain);

  Initialize(argc, argv, true, MPI_COMM_WORLD, PROB::inputs);

  Real timer_tot = amrex::second();
  Real timer_init = 0.;
  Real timer_advance = 0.;

  int max_step = -1;
  Real start_time = Real(0.0);
  Real stop_time = Real(-1.0);

  // input checks ------------------------------------------------------------
  ParmParse pp;

  pp.query("max_step", max_step);
  pp.query("stop_time", stop_time);
  if (pp.query("cfl", CNS::cfl)) {
    CNS::dt_dynamic = true;
  }
  if (pp.query("time_step", CNS::dt_constant)) {
    if (CNS::dt_dynamic) {
      amrex::Abort(
          "Simulation run parameters over-specified. Please only specify "
          "time_step or cfl");
    };
    stop_time = CNS::dt_constant * max_step;
    CNS::dt_dynamic = false;
  }

  if (start_time < Real(0.0)) {
    amrex::Abort("MUST SPECIFY a non-negative start_time");
  }

  if (max_step <= 0 || stop_time <= Real(0.0)) {
    amrex::Abort(
        "Exiting because either max_step and/or stop_time is less than or "
        "equal to 0.");
  }

  // Setup -------------------------------------------------------------------
  {
    Real timer_init = amrex::second();
    Amr amr(getLevelBld());
#ifdef AMREX_USE_GPIBM
    IBM::ib.init(
        &amr,
        NGHOST);  // this needs to happen after amr instance and before amr.init
#endif
    amr.init(start_time, stop_time);

    timer_init = amrex::second() - timer_init;
    // -------------------------------------------------------------------------

    // Time advance ------------------------------------------------------------
    amrex::Print() << " --------------------- Time advance ---------------- \n";

    timer_advance = amrex::second();
    while (amr.okToContinue() &&
               (amr.levelSteps(0) < max_step || max_step < 0) &&
               (amr.cumTime() < stop_time) ||
           stop_time < Real(0.0)) {
      //
      // Do a coarse timestep. Recursively calls timeStep()
      //
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
  // -------------------------------------------------------------------------

  timer_tot = amrex::second() - timer_tot;

  ParallelDescriptor::ReduceRealMax<Real>(
      {timer_tot, timer_init, timer_advance},
      ParallelDescriptor::IOProcessorNumber());

  amrex::Print() << "Run Time total        = " << timer_tot << "\n"
                 << "Run Time init         = " << timer_init << "\n"
                 << "Run Time advance      = " << timer_advance << "\n";

  BL_PROFILE_VAR_STOP(pmain);

  amrex::Finalize();

  return 0;
}
