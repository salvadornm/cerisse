#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Amr.H>

#include <CNS.H>

#include <IBM.H>

using namespace amrex;

amrex::LevelBld* getLevelBld ();

int main (int argc, char* argv[]) {
    amrex::Initialize(argc,argv);

    BL_PROFILE_VAR("main()", pmain);

    double timer_tot = amrex::second();
    double timer_init = 0.;
    double timer_advance = 0.;


    // Some key parameters -----------------------------------------------------
    int max_level=0;
    int  max_step = -1;
    Real strt_time = Real( 0.0);
    Real stop_time = Real(-1.0);
    {
        ParmParse pp;
        pp.query("max_level",max_level);
        pp.query("max_step",max_step);
        pp.query("strt_time",strt_time);
        pp.query("time_step",CNS::dt_glob);
        stop_time = max_step*CNS::dt_glob;
    }

    if (strt_time < Real(0.0)) {
        amrex::Abort("MUST SPECIFY a non-negative strt_time");}

    if (max_step <= 0 || stop_time <= Real(0.0)) {
        amrex::Abort("Exiting because either max_step and/or stop_time is less than or equal to 0.");}

    // Read input and setup ----------------------------------------------------
    {
        double timer_init = amrex::second();
        Amr amr(getLevelBld());
        amr.init(strt_time,stop_time);

        IBM::IB ib(amr.boxArray(),amr.DistributionMap(),2,2,max_level);
        ib.compute_markers();

        timer_init = amrex::second() - timer_init;
        exit(0);
    // -------------------------------------------------------------------------


    // Time advance ------------------------------------------------------------
    amrex::Print() << " --------------------- Time advance ---------------- \n";

        timer_advance = amrex::second();
        while ( amr.okToContinue() &&
                 (amr.levelSteps(0) < max_step || max_step < 0) &&
               (amr.cumTime() < stop_time) || stop_time < Real(0.0))

        {
            //
            // Do a coarse timestep.  Recursively calls timeStep()
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
    
    ParallelDescriptor::ReduceRealMax<double>({timer_tot, timer_init, timer_advance}, ParallelDescriptor::IOProcessorNumber());

    amrex::Print() << "Run Time total        = " << timer_tot     << "\n"
                   << "Run Time init         = " << timer_init    << "\n"
                   << "Run Time advance      = " << timer_advance << "\n";

    BL_PROFILE_VAR_STOP(pmain);

    amrex::Finalize();

    return 0;
}
