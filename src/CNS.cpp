#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParmParse.H>
#include <CNS.h>
#include <CNS_K.h>
#include <prob.h>

using namespace amrex;

bool CNS::verbose = true;
bool CNS::dt_dynamic = false;
bool CNS::ib_move = false;
bool CNS::plot_surf = false;
int CNS::nstep_screen_output = 10;
int CNS::order_rk = 2;
int CNS::stages_rk = 2;
int CNS::do_reflux = 1;
int CNS::refine_max_dengrad_lev = -1;
Real CNS::cfl = 0.0_rt;
Real CNS::dt_constant = 0.0_rt;
Real CNS::refine_dengrad = 1.0e10;

PROB::ProbClosures *CNS::h_prob_closures = nullptr;
PROB::ProbClosures *CNS::d_prob_closures = nullptr;
PROB::ProbParm *CNS::h_prob_parm = nullptr;
PROB::ProbParm *CNS::d_prob_parm = nullptr;
BCRec *CNS::h_phys_bc = nullptr;
BCRec *CNS::d_phys_bc = nullptr;

// needed for CNSBld - derived from LevelBld (abstract class, pure virtual
// functions must be implemented)

CNS::CNS() {
  // prob_rhs.init();
}

CNS::CNS(Amr &papa, int lev, const Geometry &level_geom, const BoxArray &bl,
         const DistributionMapping &dm, Real time)
    : AmrLevel(papa, lev, level_geom, bl, dm, time) {
  if (do_reflux && level > 0) {
    flux_reg.reset(new FluxRegister(grids, dmap, crse_ratio, level,PROB::ProbClosures::NCONS));
  }

#ifdef AMREX_USE_GPIBM
  IBM::ib.build_mf(grids, dmap, level);
#endif

  buildMetrics();

  // prob_rhs.init();
};

CNS::~CNS() {}
// -----------------------------------------------------------------------------

// init ------------------------------------------------------------------------

void CNS::read_params() {
  ParmParse pp("cns");

  pp.query("nstep_screen_output", nstep_screen_output);
  pp.query("verbose", verbose);

  Vector<int> lo_bc(AMREX_SPACEDIM), hi_bc(AMREX_SPACEDIM);
  pp.getarr("lo_bc", lo_bc, 0, AMREX_SPACEDIM);
  pp.getarr("hi_bc", hi_bc, 0, AMREX_SPACEDIM);
  for (int i = 0; i < AMREX_SPACEDIM; ++i) {
    h_phys_bc->setLo(i, lo_bc[i]);
    h_phys_bc->setHi(i, hi_bc[i]);
  }

  pp.query("do_reflux", do_reflux);

  if (!pp.query("order_rk", order_rk)) {
    amrex::Abort(
        "Need to specify SSPRK scheme order of accuracy, order_rk={-2, 1, 2, "
        "3}");
  }

  if (!pp.query("stages_rk", stages_rk)) {
    amrex::Abort("Need to specify SSPRK number of stages, stages_rk");
  } else {
    if (order_rk == 1 && stages_rk != 1) {
      amrex::Abort("Forward Euler number of stages must be 1");
    }

    if (order_rk == 2 && stages_rk < order_rk) {
      amrex::Abort(
          "SSPRK2 number of stages must equal or greater than order of "
          "accuracy");
    }
    if (order_rk == 3 && !(stages_rk == 4 || stages_rk == 3)) {
      amrex::Abort("SSPRK3 number of stages must equal 3 or 4");
    }
  }

#if AMREX_USE_GPIBM
  ParmParse ppib("ib");
  if (!ppib.query("move", ib_move)) {
    amrex::Abort("ib.move not specified (0=false, 1=true)");
  }
  if (!ppib.query("plot_surf", plot_surf)) {
    amrex::Abort("ib.plot_surf not specified (0=false, 1=true)");
  }
#endif

#if AMREX_USE_GPU
  amrex::Gpu::htod_memcpy(d_prob_closures, h_prob_closures,
                          sizeof(PROB::ProbClosures));
  amrex::Gpu::htod_memcpy(d_prob_parm, h_prob_parm, sizeof(PROB::ProbParm));
  amrex::Gpu::htod_memcpy(d_phys_bc, h_phys_bc, sizeof(BCRec));
#endif
}

void CNS::init(AmrLevel &old) {
  auto &oldlev = dynamic_cast<CNS &>(old);

  Real dt_new = parent->dtLevel(level);
  Real cur_time = oldlev.state[State_Type].curTime();
  Real prev_time = oldlev.state[State_Type].prevTime();
  Real dt_old = cur_time - prev_time;
  setTimeLevel(cur_time, dt_old, dt_new);

  MultiFab &S_new = get_new_data(State_Type);
  FillPatch(old, S_new, 0, cur_time, State_Type, 0,PROB::ProbClosures::NCONS);
}

void CNS::init() {
  Real dt = parent->dtLevel(level);
  Real cur_time = getLevel(level - 1).state[State_Type].curTime();
  Real prev_time = getLevel(level - 1).state[State_Type].prevTime();
  Real dt_old = (cur_time - prev_time) /
                static_cast<Real>(parent->MaxRefRatio(level - 1));
  setTimeLevel(cur_time, dt_old, dt);

  MultiFab &S_new = get_new_data(State_Type);
  FillCoarsePatch(S_new, 0, cur_time, State_Type, 0,PROB::ProbClosures::NCONS);
};

void CNS::initData() {
  BL_PROFILE("CNS::initData()");

  const auto geomdata = geom.data();
  MultiFab &S_new = get_new_data(State_Type);

  PROB::ProbClosures const *lclosures = d_prob_closures;
  PROB::ProbParm const *lprobparm = d_prob_parm;

  auto const &sma = S_new.arrays();
  amrex::ParallelFor(
      S_new, [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k) noexcept {
        prob_initdata(i, j, k, sma[box_no], geomdata, *lclosures, *lprobparm);
      });

  // TODO: Could compute primitive variables here
}

void CNS::buildMetrics() {
  // print mesh sizes
  const Real *dx = geom.CellSize();
  amrex::Print() << "Mesh size (dx,dy,dz) = ";
  amrex::Print() << AMREX_D_TERM(dx[0], << "  " << dx[1], << "  " << dx[2]) << "  \n";
}

void CNS::post_init(Real stop_time) {
  if (level > 0) {
    return;
  };

  for (int k = parent->finestLevel() - 1; k >= 0; --k) {
    getLevel(k).avgDown();
  }

  if (verbose) {
    printTotal();
  }
}
// -----------------------------------------------------------------------------

// Time-stepping ---------------------------------------------------------------
// void CNS::computeTemp(MultiFab &State, int ng) {
//   BL_PROFILE("CNS::computeTemp()");

//   PROB::ProbClosures const *lclosures = d_prob_closures;

//   // This will reset Eint and compute Temperature
// #ifdef AMREX_USE_OMP
// #pragma omp parallel if (Gpu::notInLaunchRegion())
// #endif
//   for (MFIter mfi(State, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
//     const Box &bx = mfi.growntilebox(ng);
//     auto const &sfab = State.array(mfi);

//     amrex::Abort("ComputeTemp function not written");

//     amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
//       cns_compute_temperature(i, j, k, sfab, *lclosures);
//     });
//   }
// }
// Called on level 0 only, after field data is initialised but before first
// step.
void CNS::computeInitialDt(int finest_level, int sub_cycle,
                           Vector<int> &n_cycle,  // no. of subcycling steps
                           const Vector<IntVect> &ref_ratio,
                           Vector<Real> &dt_level, Real stop_time) {
  Real dt0 = std::numeric_limits<Real>::max();
  int nfactor = 1;

  // dt at base level
  if (!dt_dynamic) {
    dt0 = dt_constant;
  } else {
    // estimate dt per level
    for (int i = 0; i <= finest_level; i++) {
      dt_level[i] = getLevel(i).estTimeStep();
    }
    // find min dt across all levels
    nfactor = 1;
    for (int i = 0; i <= finest_level; i++) {
      nfactor *= n_cycle[i];
      dt0 = std::min(dt0, nfactor * dt_level[i]);
    }
  }

  // set dt at all levels
  nfactor = 1;
  for (int i = 0; i <= finest_level; i++) {
    nfactor *= n_cycle[i];
    dt_level[i] = dt0 / nfactor;
  }

  Print() << "[computeInitialDt] Level 0 dt = " << dt0 << std::endl;
}

// Called at the end of a coarse grid timecycle or after regrid, to compute the
// dt (time step) for all levels, for the next step.
void CNS::computeNewDt(int finest_level, int /*sub_cycle*/,
                       Vector<int> &n_cycle,
                       const Vector<IntVect> & /*ref_ratio*/,
                       Vector<Real> &dt_min, Vector<Real> &dt_level,
                       Real stop_time, int post_regrid_flag) {
  // dt is constant, dt at level 0 is from inputs and dt at higher levels is
  // computed from the number of subcycles
  // ////////////////////////////////////////
  if (!dt_dynamic) {
    int nfactor = 1;
    for (int i = 0; i <= finest_level; i++) {
      nfactor *= n_cycle[i];
      dt_level[i] = dt_constant / nfactor;
    }
    return;
  }

  // if dt is dynamic //////////////////////////////////////////////////////////
  for (int i = 0; i <= finest_level; i++) {
    dt_min[i] = getLevel(i).estTimeStep();
  }

  // Limit dt
  if (post_regrid_flag == 1) {
    // Limit dt's by pre-regrid dt
    for (int i = 0; i <= finest_level; i++) {
      dt_min[i] = std::min(dt_min[i], dt_level[i]);
    }
  }
  // Limit dt's by change_max * old dt
  else {
    static Real change_max = 1.1;
    for (int i = 0; i <= finest_level; i++) {
      dt_min[i] = std::min(dt_min[i], change_max * dt_level[i]);
    }
  }

  // Find the minimum over all levels
  Real dt0 = std::numeric_limits<Real>::max();
  int nfactor = 1;
  for (int i = 0; i <= finest_level; i++) {
    nfactor *= n_cycle[i];
    dt0 = std::min(dt0, nfactor * dt_min[i]);
  }

  // Limit dt0 by the value of stop_time.
  const Real eps = 0.001_rt * dt0;
  Real cur_time = state[State_Type].curTime();
  if (stop_time >= 0.0_rt) {
    if ((cur_time + dt0) > (stop_time - eps)) {
      dt0 = stop_time - cur_time;
    }
  }

  // Set dt at all levels
  nfactor = 1;
  for (int i = 0; i <= finest_level; i++) {
    nfactor *= n_cycle[i];
    dt_level[i] = dt0 / nfactor;
  }

  Print() << "[computeNewDt] Level 0 dt = " << dt0 << std::endl;
}

Real CNS::estTimeStep() {
  BL_PROFILE("CNS::estTimeStep()");

  const auto dx = geom.CellSizeArray();
  const MultiFab &S = get_new_data(State_Type);
  PROB::ProbClosures const *lclosures = d_prob_closures;

  Real estdt = amrex::ReduceMin(
      S, 0,
      [=] AMREX_GPU_HOST_DEVICE(Box const &bx, Array4<Real const> const &fab)
          -> Real { return cns_estdt(bx, fab, dx, *lclosures); });

  estdt *= cfl;
  ParallelDescriptor::ReduceRealMin(estdt);

  return estdt;
}

void CNS::post_timestep(int /* iteration*/) {
  BL_PROFILE("post_timestep");

  if (do_reflux && level < parent->finestLevel()) {
    MultiFab &S = get_new_data(State_Type);
    CNS &fine_level = getLevel(level + 1);
    fine_level.flux_reg->Reflux(S, Real(1.0), 0, 0, PROB::ProbClosures::NCONS, geom);
  }

  if (level < parent->finestLevel()) {
    avgDown();
  }

  if (verbose && ((this->nStep() % nstep_screen_output) == 0)) {
    printTotal();
  }
}

void CNS::postCoarseTimeStep(Real time) {
#if AMREX_USE_GPIBM
  // if (ib_move) {
  //   IBM::ib.moveGeom();
  //   // reallocate variables?
  //   // Print() << parent->finestLevel() << std::endl;
  //   for (int lev=0; lev <= parent->finestLevel(); lev++) {
  //     IBM::ib.computeMarkers(0);
  //     IBM::ib.initialiseGPs(0);
  //   }
  // }
#endif
}
// -----------------------------------------------------------------------------

// Gridding -------------------------------------------------------------------
// Called for each level from 0,1...nlevs-1
void CNS::post_regrid(int lbase, int new_finest) {
#ifdef AMREX_USE_GPIBM
  IBM::ib.destroy_mf(level);
  IBM::ib.build_mf(grids, dmap, level);
  IBM::ib.computeMarkers(level);
  IBM::ib.initialiseGPs(level);
#endif
}

void CNS::errorEst(TagBoxArray &tags, int /*clearval*/, int /*tagval*/,
                   Real time, int /*n_error_buf*/, int /*ngrow*/) {
  // MF without ghost points filled (why?)
  MultiFab sdata(get_new_data(State_Type).boxArray(),
                 get_new_data(State_Type).DistributionMap(), PROB::ProbClosures::NCONS, PROB::ProbClosures::NGHOST,
                 MFInfo(), Factory());

  // filling ghost points (copied from PeleC)
  const Real cur_time = state[State_Type].curTime();
  FillPatch(*this, sdata, PROB::ProbClosures::NGHOST, cur_time, State_Type, 0, PROB::ProbClosures::NCONS);
  const auto geomdata = geom.data();

#ifdef AMREX_USE_GPIBM
  // call function from cns_prob
  auto &ibdata = (*IBM::ib.bmf_a[level]);
#endif
  for (MFIter mfi(tags, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
    const Box &bx = mfi.tilebox();
    auto const &tagfab = tags.array(mfi);
    auto const &sdatafab = sdata.array(mfi);
#ifdef AMREX_USE_GPIBM
    auto const &ibfab = ibdata.array(mfi);
#endif
    int lev = level;
    int nt_lev = nStep();
    PROB::ProbParm const *lprobparm = d_prob_parm;

    ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
#ifdef AMREX_USE_GPIBM
      // call function from cns_prob
      user_tagging(i, j, k, nt_lev, tagfab, sdatafab, ibfab, geomdata,
                   *lprobparm, lev);
#else
    user_tagging(i, j, k, nt_lev, tagfab, sdatafab, geomdata ,*lprobparm, lev);
#endif
    });
  }
}

// TODO: Add restarts
// -----------------------------------------------------------------------------

void CNS::avgDown() {
  BL_PROFILE("CNS::avgDown()");

  if (level == parent->finestLevel()) return;

  auto &fine_lev = getLevel(level + 1);

  MultiFab &S_crse = get_new_data(State_Type);
  MultiFab &S_fine = fine_lev.get_new_data(State_Type);

  amrex::average_down(S_fine, S_crse, fine_lev.geom, geom, 0, S_fine.nComp(),
                      parent->refRatio(level));
}

void CNS::printTotal() const {
  // Get conservatives multifab
  const MultiFab& consmf = get_new_data(State_Type);

  // Compute peicewise constant sum of conserved variables
  std::array<Real, PROB::ProbClosures::NCONS> tot;
  for (int comp = 0; comp < PROB::ProbClosures::NCONS; ++comp) {
    tot[comp] = consmf.sum(comp, true) * geom.ProbSize();
  }

  // Communicate across processors
  ParallelDescriptor::ReduceRealSum(tot.data(), PROB::ProbClosures::NCONS,
                                    ParallelDescriptor::IOProcessorNumber());
  // Print
  Vector<std::string> names= PROB::ProbClosures::get_cons_vars_names();
  for (int comp = 0; comp < PROB::ProbClosures::NCONS; ++comp) {
    amrex::Print().SetPrecision(17) << "   Total " << names[comp] << " = " << tot[comp] << "\n";
  }
}

void CNS::variableCleanUp() {
  delete h_prob_closures;
  delete h_phys_bc;

#ifdef AMREX_USE_GPU
  The_Arena()->free(d_prob_closures);
  The_Arena()->free(d_phys_bc);
#endif
  desc_lst.clear();
  derive_lst.clear();
}

// Plotting
//------------------------------------------------------------------------------
void CNS::writePlotFile(const std::string &dir, std::ostream &os,
                        VisMF::How how) {
  int i, n;
  //
  // The list of indices of State to write to plotfile.
  // first component of pair is state_type,
  // second component of pair is component # within the state_type
  //
  std::vector<std::pair<int, int>> plot_var_map;
  for (int typ = 0; typ < desc_lst.size(); typ++) {
    for (int comp = 0; comp < desc_lst[typ].nComp(); comp++) {
      if (parent->isStatePlotVar(desc_lst[typ].name(comp)) &&
          desc_lst[typ].getType() == IndexType::TheCellType()) {
        plot_var_map.push_back(std::pair<int, int>(typ, comp));
      }
    }
  }

  int num_derive = 0;
  std::vector<std::string> derive_names;
  const std::list<DeriveRec> &dlist = derive_lst.dlist();
  for (auto const &d : dlist) {
    if (parent->isDerivePlotVar(d.name())) {
      derive_names.push_back(d.name());
      num_derive += d.numDerive();
    }
  }

  int n_data_items = plot_var_map.size() + num_derive;

//----------------------------------------------------------------------modified
#ifdef AMREX_USE_GPIBM
  n_data_items += 2;
#endif
  //------------------------------------------------------------------------------

  // get the time from the first State_Type
  // if the State_Type is ::Interval, this will get t^{n+1/2} instead of t^n
  Real cur_time = state[0].curTime();

  if (level == 0 && ParallelDescriptor::IOProcessor()) {
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
    for (i = 0; i < static_cast<int>(plot_var_map.size()); i++) {
      int typ = plot_var_map[i].first;
      int comp = plot_var_map[i].second;
      os << desc_lst[typ].name(comp) << '\n';
    }

    // derived
    for (auto const &dname : derive_names) {
      const DeriveRec *rec = derive_lst.get(dname);
      for (i = 0; i < rec->numDerive(); ++i) {
        os << rec->variableName(i) << '\n';
      }
    }

//----------------------------------------------------------------------modified
#ifdef AMREX_USE_GPIBM
    os << "sld\n";
    os << "ghs\n";
#endif
    //------------------------------------------------------------------------------

    os << AMREX_SPACEDIM << '\n';
    os << parent->cumTime() << '\n';
    int f_lev = parent->finestLevel();
    os << f_lev << '\n';
    for (i = 0; i < AMREX_SPACEDIM; i++) os << Geom().ProbLo(i) << ' ';
    os << '\n';
    for (i = 0; i < AMREX_SPACEDIM; i++) os << Geom().ProbHi(i) << ' ';
    os << '\n';
    for (i = 0; i < f_lev; i++) os << parent->refRatio(i)[0] << ' ';
    os << '\n';
    for (i = 0; i <= f_lev; i++) os << parent->Geom(i).Domain() << ' ';
    os << '\n';
    for (i = 0; i <= f_lev; i++) os << parent->levelSteps(i) << ' ';
    os << '\n';
    for (i = 0; i <= f_lev; i++) {
      for (int k = 0; k < AMREX_SPACEDIM; k++)
        os << parent->Geom(i).CellSize()[k] << ' ';
      os << '\n';
    }
    os << (int)Geom().Coord() << '\n';
    os << "0\n";  // Write bndry data.
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
  if (!FullPath.empty() && FullPath[FullPath.size() - 1] != '/') {
    FullPath += '/';
  }
  FullPath += sLevel;
  //
  // Only the I/O processor makes the directory if it doesn't already exist.
  //
  if (!levelDirectoryCreated) {
    if (ParallelDescriptor::IOProcessor()) {
      if (!amrex::UtilCreateDirectory(FullPath, 0755)) {
        amrex::CreateDirectoryFailed(FullPath);
      }
    }
    // Force other processors to wait until directory is built.
    ParallelDescriptor::Barrier();
  }

  if (ParallelDescriptor::IOProcessor()) {
    os << level << ' ' << grids.size() << ' ' << cur_time << '\n';
    os << parent->levelSteps(level) << '\n';

    for (i = 0; i < grids.size(); ++i) {
      RealBox gridloc = RealBox(grids[i], geom.CellSize(), geom.ProbLo());
      for (n = 0; n < AMREX_SPACEDIM; n++)
        os << gridloc.lo(n) << ' ' << gridloc.hi(n) << '\n';
    }
    //
    // The full relative pathname of the MultiFabs at this level.
    // The name is relative to the Header file containing this name.
    // It's the name that gets written into the Header.
    //
    if (n_data_items > 0) {
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
  int cnt = 0;
  const int nGrow = 0;
  MultiFab plotMF(grids, dmap, n_data_items, nGrow, MFInfo(), Factory());
  MultiFab *this_dat = 0;
  //
  // Cull data from state variables -- use no ghost cells.
  //
  for (i = 0; i < static_cast<int>(plot_var_map.size()); i++) {
    int typ = plot_var_map[i].first;
    int comp = plot_var_map[i].second;
    this_dat = &state[typ].newData();
    MultiFab::Copy(plotMF, *this_dat, comp, cnt, 1, nGrow);
    cnt++;
  }

  // derived
  if (derive_names.size() > 0) {
    for (auto const &dname : derive_names) {
      derive(dname, cur_time, plotMF, cnt);
      cnt += derive_lst.get(dname)->numDerive();
    }
  }

//----------------------------------------------------------------------modified
#ifdef AMREX_USE_GPIBM
  plotMF.setVal(0.0_rt, cnt, 2, nGrow);
  IBM::ib.bmf_a[level]->copytoRealMF(plotMF, 0, cnt);
#endif
  //------------------------------------------------------------------------------

  //
  // Use the Full pathname when naming the MultiFab.
  //
  std::string TheFullPath = FullPath;
  TheFullPath += BaseName;
  if (AsyncOut::UseAsyncOut()) {
    VisMF::AsyncWrite(plotMF, TheFullPath);
  } else {
    VisMF::Write(plotMF, TheFullPath, how, true);
  }

  levelDirectoryCreated = false;  // ---- now that the plotfile is finished
}

// This is called once per level on write timestep.
void CNS::writePlotFilePost(const std::string &dir, std::ostream &os) {
#if AMREX_USE_GPIBM
  // write geometry -- each proc holds same geom, even with ib_move, so no
  // communication is required. if (plot_surf) { Print() << "Writing surface
  // data" << std::endl; if (ioproc=0) ib.writeGeom()

  // IBM::ib.computeSurf(this->level); // computed at each level. From low to
  // high. Print() << "Computed surface data" << std::endl;

  // Print() << "Writing surface data" << std::endl;
  // if (ioproc=0) ib.writeSurf()
  // }
#endif
}