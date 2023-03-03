#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParmParse.H>
#include <CNS.H>
#include <CNS_K.H>
#include <cns_prob.H>
#include <CNS_parm.H>
#include <IBM.H>

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
  Real time, int /*n_error_buf*/, int /*ngrow*/) {

  // MF without ghost points filled (why?)
  MultiFab data( get_new_data(State_Type).boxArray(), get_new_data(State_Type).DistributionMap(), NUM_STATE, 1, MFInfo(), Factory());

  // filling ghost points (copied from PeleC)
  const Real cur_time = state[State_Type].curTime();
  FillPatch(*this, data, data.nGrow(), cur_time, State_Type, 0, NUM_STATE, 0);

  // IBM::IBMultiFab *ibmf = IBM::ib.mfa->at(level); TODO - pass ibmf reference to tagging

  // call function from cns_prob
  tagging(tags, data, level);

  // -------------------------------- Monal 03/03/23
  // amrex::Print() << state[0].descriptor()->name(0) << std::endl;
  //
  // In tagging(...) retrieving data MultiFabs seem to have ghost points not filled - atleast straight after initialisation.
  // I have tried the following:
  // MultiFab &data=get_new_data(0); 
  // MultiFab& data = state.at(0).newData();
  // both give erroneous computation error.
  //
  // PeleC creates a local copy of the MultiFab rather than passing through the master data with:
  // MultiFab data( get_new_data(State_Type).boxArray(), get_new_data(State_Type).DistributionMap(), NUM_STATE, 1, MFInfo(), Factory());
  // However, this is not necessary.
  //---------------------------------
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

// Plotting
//------------------------------------------------------------------------------
void CNS::writePlotFile (const std::string& dir,
                         std::ostream&      os,
                         VisMF::How         how)
{
    // amrex::Print() << "custom plot " << " (level" << level << ") "<< dir << std::endl;
    // amrex::Print() << parent->levelSteps(level) << std::endl;
    int i, n;
    //
    // The list of indices of State to write to plotfile.
    // first component of pair is state_type,
    // second component of pair is component # within the state_type
    //
    std::vector<std::pair<int,int> > plot_var_map;
    for (int typ = 0; typ < desc_lst.size(); typ++)
    {
        for (int comp = 0; comp < desc_lst[typ].nComp();comp++)
        {
            if (parent->isStatePlotVar(desc_lst[typ].name(comp)) &&
                desc_lst[typ].getType() == IndexType::TheCellType())
            {
                plot_var_map.push_back(std::pair<int,int>(typ,comp));
            }
        }
    }

    int num_derive = 0;
    std::vector<std::string> derive_names;
    const std::list<DeriveRec>& dlist = derive_lst.dlist();
    for (auto const& d : dlist)
    {
        if (parent->isDerivePlotVar(d.name()))
        {
            derive_names.push_back(d.name());
            num_derive += d.numDerive();
        }
    }

    int n_data_items = plot_var_map.size() + num_derive;

//----------------------------------------------------------------------modified
// #ifdef AMREX_USE_EB
    // if (EB2::TopIndexSpaceIfPresent()) {
        n_data_items += 1;
    // }
// #endif
//------------------------------------------------------------------------------

    // get the time from the first State_Type
    // if the State_Type is ::Interval, this will get t^{n+1/2} instead of t^n
    Real cur_time = state[0].curTime();

    if (level == 0 && ParallelDescriptor::IOProcessor())
    {
        //
        // The first thing we write out is the plotfile type.
        //
        os << thePlotFileType() << '\n';

        if (n_data_items == 0)
            amrex::Error("Must specify at least one valid data item to plot");

        os << n_data_items << '\n';

        //
        // Names of variables
        //
        for (i =0; i < static_cast<int>(plot_var_map.size()); i++)
        {
            int typ = plot_var_map[i].first;
            int comp = plot_var_map[i].second;
            os << desc_lst[typ].name(comp) << '\n';
        }

        // derived
        for (auto const& dname : derive_names) {
            const DeriveRec* rec = derive_lst.get(dname);
            for (i = 0; i < rec->numDerive(); ++i) {
                os << rec->variableName(i) << '\n';
            }
        }

//----------------------------------------------------------------------modified
// #ifdef AMREX_USE_EB
        // if (EB2::TopIndexSpaceIfPresent()) {
            os << "sld\n";
        // }
// #endif
//------------------------------------------------------------------------------

        os << AMREX_SPACEDIM << '\n';
        os << parent->cumTime() << '\n';
        int f_lev = parent->finestLevel();
        os << f_lev << '\n';
        for (i = 0; i < AMREX_SPACEDIM; i++)
            os << Geom().ProbLo(i) << ' ';
        os << '\n';
        for (i = 0; i < AMREX_SPACEDIM; i++)
            os << Geom().ProbHi(i) << ' ';
        os << '\n';
        for (i = 0; i < f_lev; i++)
            os << parent->refRatio(i)[0] << ' ';
        os << '\n';
        for (i = 0; i <= f_lev; i++)
            os << parent->Geom(i).Domain() << ' ';
        os << '\n';
        for (i = 0; i <= f_lev; i++)
            os << parent->levelSteps(i) << ' ';
        os << '\n';
        for (i = 0; i <= f_lev; i++)
        {
            for (int k = 0; k < AMREX_SPACEDIM; k++)
                os << parent->Geom(i).CellSize()[k] << ' ';
            os << '\n';
        }
        os << (int) Geom().Coord() << '\n';
        os << "0\n"; // Write bndry data.

    }
    // Build the directory to hold the MultiFab at this level.
    // The name is relative to the directory containing the Header file.
    //
    static const std::string BaseName = "/Cell";
    char buf[64];
    snprintf(buf, sizeof buf, "Level_%d", level);
    std::string sLevel = buf;
    //
    // Now for the full pathname of that directory.
    //
    std::string FullPath = dir;
    if ( ! FullPath.empty() && FullPath[FullPath.size()-1] != '/')
    {
        FullPath += '/';
    }
    FullPath += sLevel;
    //
    // Only the I/O processor makes the directory if it doesn't already exist.
    //
    if ( ! levelDirectoryCreated) {
        if (ParallelDescriptor::IOProcessor()) {
            if ( ! amrex::UtilCreateDirectory(FullPath, 0755)) {
                amrex::CreateDirectoryFailed(FullPath);
            }
        }
        // Force other processors to wait until directory is built.
        ParallelDescriptor::Barrier();
    }

    if (ParallelDescriptor::IOProcessor())
    {
        os << level << ' ' << grids.size() << ' ' << cur_time << '\n';
        os << parent->levelSteps(level) << '\n';

        for (i = 0; i < grids.size(); ++i)
        {
            RealBox gridloc = RealBox(grids[i],geom.CellSize(),geom.ProbLo());
            for (n = 0; n < AMREX_SPACEDIM; n++)
                os << gridloc.lo(n) << ' ' << gridloc.hi(n) << '\n';
        }
        //
        // The full relative pathname of the MultiFabs at this level.
        // The name is relative to the Header file containing this name.
        // It's the name that gets written into the Header.
        //
        if (n_data_items > 0)
        {
            std::string PathNameInHeader = sLevel;
            PathNameInHeader += BaseName;
            os << PathNameInHeader << '\n';
        }

//----------------------------------------------------------------------modified
// #ifdef AMREX_USE_EB
        // if (EB2::TopIndexSpaceIfPresent()) {
        //     volfrac threshold for amrvis
        //     if (level == parent->finestLevel()) {
        //         for (int lev = 0; lev <= parent->finestLevel(); ++lev) {
        //             os << "1.0e-6\n";
        //         }
        //     }
        // }
// #endif
//------------------------------------------------------------------------------
    }
    //
    // We combine all of the multifabs -- state, derived, etc -- into one
    // multifab -- plotMF.
    int       cnt   = 0;
    const int nGrow = 0;
    MultiFab  plotMF(grids,dmap,n_data_items,nGrow,MFInfo(),Factory());
    MultiFab* this_dat = 0;
    //
    // Cull data from state variables -- use no ghost cells.
    //
    for (i = 0; i < static_cast<int>(plot_var_map.size()); i++)
    {
        int typ  = plot_var_map[i].first;
        int comp = plot_var_map[i].second;
        this_dat = &state[typ].newData();
        MultiFab::Copy(plotMF,*this_dat,comp,cnt,1,nGrow);
        cnt++;
    }

    // derived
    if (derive_names.size() > 0)
    {
        for (auto const& dname : derive_names)
        {
            derive(dname, cur_time, plotMF, cnt);
            cnt += derive_lst.get(dname)->numDerive();
        }
    }

//----------------------------------------------------------------------modified
// #ifdef AMREX_USE_EB
    // if (EB2::TopIndexSpaceIfPresent()) {
        plotMF.setVal(0.0, cnt, 1, nGrow);
        // auto factory = static_cast<EBFArrayBoxFactory*>(m_factory.get());

if (parent->levelSteps(level)>0) {
        IBM::IBMultiFab* ibmf = IBM::ib.mfa->at(level);
        ibmf->copytoRealMF(plotMF,0,cnt);}
        // MultiFab::Copy(plotMF,IB,0,cnt,1,nGrow);
    // }
// #endif
//------------------------------------------------------------------------------

    //
    // Use the Full pathname when naming the MultiFab.
    //
    std::string TheFullPath = FullPath;
    TheFullPath += BaseName;
    if (AsyncOut::UseAsyncOut()) {
        VisMF::AsyncWrite(plotMF,TheFullPath);
    } else {
        VisMF::Write(plotMF,TheFullPath,how,true);
    }

    levelDirectoryCreated = false;  // ---- now that the plotfile is finished
}
