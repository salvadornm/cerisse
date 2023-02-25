#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParmParse.H>
#include <CNS.H>
#include <CNS_K.H>
#include <cns_prob.H>
#include <CNS_parm.H>

using namespace amrex;

constexpr int CNS::NUM_GROW;
BCRec     CNS::phys_bc;

int       CNS::verbose = 0;
IntVect   CNS::hydro_tile_size {AMREX_D_DECL(1024,16,16)};
Real      CNS::cfl       = 0.3;
Real      CNS::dt_glob = 0.0;
int       CNS::do_reflux = 1;
int       CNS::refine_max_dengrad_lev   = -1;
Real      CNS::refine_dengrad           = 1.0e10;

Real      CNS::gravity = 0.0;

// needed for CNSBld - derived from LevelBld (abstract class, pure virtual functions must be implemented) 

CNS::CNS () {}

CNS::CNS (Amr&            papa,
          int             lev,
          const Geometry& level_geom,
          const BoxArray& bl,
          const DistributionMapping& dm,
          Real            time)
    : AmrLevel(papa,lev,level_geom,bl,dm,time)
{
    if (do_reflux && level > 0) {
        flux_reg.reset(new FluxRegister(grids,dmap,crse_ratio,level,NUM_STATE));
    }

    buildMetrics();
}

CNS::~CNS () {}
// -----------------------------------------------------------------------------

// init ------------------------------------------------------------------------

void CNS::read_params () {
    ParmParse pp("cns");

    // pp.query("v", verbose);

    Vector<int> tilesize(AMREX_SPACEDIM);
    if (pp.queryarr("hydro_tile_size", tilesize, 0, AMREX_SPACEDIM))
    {
        for (int i=0; i<AMREX_SPACEDIM; i++) hydro_tile_size[i] = tilesize[i];
    }

    pp.query("cfl", cfl);
    pp.query("time_step",dt_glob);
    // pp.query("n_cycle",n_cycle);

    Vector<int> lo_bc(AMREX_SPACEDIM), hi_bc(AMREX_SPACEDIM);
    pp.getarr("lo_bc", lo_bc, 0, AMREX_SPACEDIM);
    pp.getarr("hi_bc", hi_bc, 0, AMREX_SPACEDIM);
    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
        phys_bc.setLo(i,lo_bc[i]);
        phys_bc.setHi(i,hi_bc[i]);
    }

    pp.query("do_reflux", do_reflux);

    pp.query("refine_max_dengrad_lev", refine_max_dengrad_lev);
    pp.query("refine_dengrad", refine_dengrad);

    pp.query("gravity", gravity);

    // pp.query("eos_gamma", h_parm->eos_gamma);

    parm->Initialize();
}

void CNS::init (AmrLevel& old) {
    auto& oldlev = dynamic_cast<CNS&>(old);

    Real dt_new    = parent->dtLevel(level);
    Real cur_time  = oldlev.state[State_Type].curTime();
    Real prev_time = oldlev.state[State_Type].prevTime();
    Real dt_old    = cur_time - prev_time;
    setTimeLevel(cur_time,dt_old,dt_new);

    MultiFab& S_new = get_new_data(State_Type);
    FillPatch(old,S_new,0,cur_time,State_Type,0,NUM_STATE);
}

void CNS::init() {
    Real dt        = parent->dtLevel(level);
    Real cur_time  = getLevel(level-1).state[State_Type].curTime();
    Real prev_time = getLevel(level-1).state[State_Type].prevTime();
    Real dt_old = (cur_time - prev_time)/static_cast<Real>(parent->MaxRefRatio(level-1));
    setTimeLevel(cur_time,dt_old,dt);

    MultiFab& S_new = get_new_data(State_Type);
    FillCoarsePatch(S_new, 0, cur_time, State_Type, 0, NUM_STATE);
};

void CNS::initData ()
{
    BL_PROFILE("CNS::initData()");

    const auto geomdata = geom.data();
    MultiFab& S_new = get_new_data(State_Type);

    // Parm const* lparm = d_parm;
    // ProbParm const* lprobparm = d_prob_parm;

    auto const& sma = S_new.arrays();
    amrex::ParallelFor(S_new,
    [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
    {
        cns_initdata(i, j, k, sma[box_no], geomdata, *parm, *prob_parm);
    });

    // Compute the initial temperature (will override what was set in initdata)
    // computeTemp(S_new,0);

    MultiFab& C_new = get_new_data(Cost_Type);
    C_new.setVal(1.0);
}

void CNS::buildMetrics ()
{
    // print mesh sizes
    const Real* dx = geom.CellSize();
    amrex::Print() << "Mesh size (dx,dy,dz) = ";
    amrex::Print() << dx[0] << "  "
                   << dx[1] << "  "
                   << dx[2] << "  \n";
}

void CNS::post_init (Real) {
    if (level > 0) return;
    for (int k = parent->finestLevel()-1; k >= 0; --k) {
        getLevel(k).avgDown();
    }

    // if (verbose >= 2) {
        printTotal();
    // }
}

// -----------------------------------------------------------------------------


// Time-stepping ---------------------------------------------------------------
void CNS::computeTemp (MultiFab& State, int ng)
{
    BL_PROFILE("CNS::computeTemp()");

    Parm const* lparm = parm;

    // This will reset Eint and compute Temperature
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(State,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.growntilebox(ng);
        auto const& sfab = State.array(mfi);

        amrex::ParallelFor(bx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            cns_compute_temperature(i,j,k,sfab,*lparm);
        });
    }
}

Real CNS::initialTimeStep () {
    return CNS::dt_glob;
}

void CNS::computeInitialDt (int                finest_level,
                       int                    /*sub_cycle*/,
                       Vector<int>&           n_cycle,
                       const Vector<IntVect>& /*ref_ratio*/,
                       Vector<Real>&          dt_level,
                       Real                   stop_time) {
    //
    // Grids have been constructed, compute dt for all levels.
    //
    // if (level > 0) {
    //     return;
    // }

    // Real dt_0 = std::numeric_limits<Real>::max();
    // int n_factor = 1;
    // for (int i = 0; i <= finest_level; i++)
    // {
    //     dt_level[i] = getLevel(i).initialTimeStep();
    //     n_factor   *= n_cycle[i];
    //     dt_0 = std::min(dt_0,n_factor*dt_level[i]);
    // }
    int n_factor = 1;
    dt_level[0] = getLevel(0).initialTimeStep();
    for (int i = 0; i <= finest_level; i++)
    {
        n_factor *= n_cycle[i];
        dt_level[i] = dt_level[0]/n_factor;
        // amrex::Print() << "level i="<< i << ", dt ="<< dt_level[i] << "\n";
    }

    //
    // Limit dt's by the value of stop_time.
    //
    // const Real eps = 0.001*dt_0;
    // Real cur_time  = state[State_Type].curTime();
    // if (stop_time >= 0.0) {
    //     if ((cur_time + dt_0) > (stop_time - eps))
    //         dt_0 = stop_time - cur_time;
    // }


}

/////////////////// SAME AS computeInitialDt for now ///////////////////////////
void CNS::computeNewDt (int                    finest_level,
                   int                    /*sub_cycle*/,
                   Vector<int>&           n_cycle,
                   const Vector<IntVect>& /*ref_ratio*/,
                   Vector<Real>&          dt_min,
                   Vector<Real>&          dt_level,
                   Real                   stop_time,
                   int                    post_regrid_flag) {
    //
    // Grids have been constructed, compute dt for all levels.
    //
    // if (level > 0) {
    //     return;
    // }

    // Real dt_0 = std::numeric_limits<Real>::max();
    // int n_factor = 1;
    // for (int i = 0; i <= finest_level; i++)
    // {
    //     dt_level[i] = getLevel(i).initialTimeStep();
    //     n_factor   *= n_cycle[i];
    //     dt_0 = std::min(dt_0,n_factor*dt_level[i]);
    // }
    int n_factor = 1;
    dt_level[0] = getLevel(0).initialTimeStep();
    for (int i = 0; i <= finest_level; i++)
    {
        n_factor *= n_cycle[i];
        dt_level[i] = dt_level[0]/n_factor;
        // amrex::Print() << "level i="<< i << ", dt ="<< dt_level[i] << "\n";
    }

}

void CNS::post_timestep (int /*iteration*/) {
    BL_PROFILE("post_timestep");

    if (do_reflux && level < parent->finestLevel()) {
        MultiFab& S = get_new_data(State_Type);
        CNS& fine_level = getLevel(level+1);
        fine_level.flux_reg->Reflux(S, Real(1.0), 0, 0, NUM_STATE, geom);
    }

    if (level < parent->finestLevel()) {
        avgDown();
    }
}
// -----------------------------------------------------------------------------

// Gridding -------------------------------------------------------------------

void CNS::post_regrid (int /*lbase*/, int /*new_finest*/) { }


void CNS::errorEst (TagBoxArray& tags, int /*clearval*/, int /*tagval*/,
  amrex::Real time, int /*n_error_buf*/, int /*ngrow*/) {

  // get state MF

  // Without Ghost points
  // const amrex::MultiFab& data = get_new_data(State_Type);

  // With Ghost points -------copied from PeleC. Can we avoid the Fillpatch?
  amrex::MultiFab data( get_new_data(State_Type).boxArray(), get_new_data(State_Type).DistributionMap(), NUM_STATE, 1, amrex::MFInfo(), Factory());

  const amrex::Real cur_time = state[State_Type].curTime();
  FillPatch(*this, data, data.nGrow(), cur_time, State_Type, Density, NUM_STATE, 0);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(tags,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        auto const& tagfab = tags.array(mfi);
        auto const& datafab = data.array(mfi);

        amrex::ParallelFor(bx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          // Temporary tagging on density in x only
          Real drhox = amrex::Math::abs(datafab(i+1,j,k,0) - datafab(i-1,j,k,0))/datafab(i,j,k,0);
          tagfab(i,j,k) = drhox > 0.5f;

          // amrex::Print() << "i,j,k         " << i << " " << j << " " << k << " "<< std::endl;
          // amrex::Print() << "tag(i,j,k)    " << int(tagfab(i,j,k)) << std::endl;
          // amrex::Print() << "data(i,j,k,0) " << datafab(i,j,k,0) << std::endl;
          // amrex::Print() << "drhox         " << drhox << std::endl;
          // amrex::Print() << "------------- " << std::endl;
        });
    }
}
// -----------------------------------------------------------------------------

// misc ------------------------------------------------------------------------

void CNS::avgDown () {
    BL_PROFILE("CNS::avgDown()");

    if (level == parent->finestLevel()) return;

    auto& fine_lev = getLevel(level+1);

    MultiFab& S_crse =          get_new_data(State_Type);
    MultiFab& S_fine = fine_lev.get_new_data(State_Type);

    amrex::average_down(S_fine, S_crse, fine_lev.geom, geom,
                        0, S_fine.nComp(), parent->refRatio(level));

    const int nghost = 0;
    /////////////////////////TODO///////////////////////////
    // computeTemp(S_crse, nghost);
}

void CNS::printTotal () const {
    const MultiFab& S_new = get_new_data(State_Type);
    std::array<Real,5> tot;
    for (int comp = 0; comp < 5; ++comp) {
        tot[comp] = S_new.sum(comp,true) * geom.ProbSize();
    }
#ifdef BL_LAZY
    Lazy::QueueReduction( [=] () mutable {
#endif
            ParallelDescriptor::ReduceRealSum(tot.data(), 5, ParallelDescriptor::IOProcessorNumber());
            amrex::Print().SetPrecision(17) << "\n[CNS] Total mass       is " << tot[0] << "\n"
                                            <<   "      Total x-momentum is " << tot[1] << "\n"
                                            <<   "      Total y-momentum is " << tot[2] << "\n"
                                            <<   "      Total z-momentum is " << tot[3] << "\n"
                                            <<   "      Total energy     is " << tot[4] << "\n";
#ifdef BL_LAZY
        });
#endif
}

