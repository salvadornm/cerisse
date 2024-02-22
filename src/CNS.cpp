#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParmParse.H>
#include <CNS.h>
#include <CNS_K.h>
#include <prob.h>

#ifdef AMREX_USE_GPIBM
#include <IBM.h>
#endif

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

CNS::CNS() {}

CNS::CNS(Amr &papa, int lev, const Geometry &level_geom, const BoxArray &bl,
         const DistributionMapping &dm, Real time)
    : AmrLevel(papa, lev, level_geom, bl, dm, time) {
  if (do_reflux && level > 0) {
    flux_reg.reset(new FluxRegister(grids, dmap, crse_ratio, level,PROB::ProbClosures::NCONS));
  }

#ifdef AMREX_USE_GPIBM
  IBM::ib.buildMFs(grids, dmap, level);
#endif

  buildMetrics();
}

CNS::~CNS() {}
// -----------------------------------------------------------------------------

// init ------------------------------------------------------------------------

void CNS::read_params() {
  ParmParse pp("cns");

  pp.query("screen_output", nstep_screen_output);
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
  amrex::Print() << dx[0] << "  " << dx[1] << "  " << dx[2] << "  \n";
}

void CNS::post_init(Real stop_time) {
  if (level > 0) {
    return;
  };

  for (int k = parent->finestLevel() - 1; k >= 0; --k) {
    getLevel(k).avgDown();
  }

  if (verbose) {
    // printTotal();
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

  if (verbose && this->nStep() % nstep_screen_output == 0) {
    // printTotal();
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
  IBM::ib.destroyMFs(level);
  IBM::ib.buildMFs(grids, dmap, level);
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
  IBM::IBMultiFab &ibdata = (*IBM::ib.ibMFa[level]);
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
//   const MultiFab &S_new = get_new_data(State_Type);
//   std::array<Real, NCONS> tot;
//   std::array<Real, NPRIM> prims_max, prims_min;

//   // NEED to put these lines in GPU launch region?
//   for (int comp = 0; comp < NCONS; ++comp) {
//     tot[comp] = S_new.sum(comp, true) * geom.ProbSize();
//   }

//   for (int comp = 0; comp < NPRIM; ++comp) {
//     prims_max[comp] = Vprimsmf[level].max(comp, 0, true);
//     prims_min[comp] = Vprimsmf[level].min(comp, 0, true);
//   }

// #ifdef BL_LAZY
//   Lazy::QueueReduction([=]() mutable {
// #endif
//     ParallelDescriptor::ReduceRealSum(tot.data(), NCONS,
//                                       ParallelDescriptor::IOProcessorNumber());

//     ParallelDescriptor::ReduceRealMax(prims_max.data(), NPRIM,
//                                       ParallelDescriptor::IOProcessorNumber());

//     ParallelDescriptor::ReduceRealMin(prims_min.data(), NPRIM,
//                                       ParallelDescriptor::IOProcessorNumber());

//     // compute convective CFL
//     const auto dx = geom.CellSizeArray();
//     const Real dt = parent->dtLevel(level);
//     const MultiFab &primsmf = Vprimsmf[level];
//     PROB::ProbClosures const &lclosures = *d_prob_closures;
//     Array2D<Real, 0, 2, 0, 2> *arrayCFL;

// #if AMREX_USE_GPU
//     arrayCFL = (Array2D<Real, 0, 2, 0, 2> *)The_Arena()->alloc(
//         sizeof(Array2D<Real, 0, 2, 0, 2>));
// #else
//   arrayCFL = new Array2D<Real, 0, 2, 0, 2>{};
// #endif
//     // We cannot modify variables defined outside of the lambda function in the
//     // lambda function (not even reals and ints). AMReX does not allow mutable
//     // keyword. This is why we need to call functions.

//     for (MFIter mfi(primsmf, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
//       auto const &prims = primsmf.array(mfi);
//       const Box &bx = mfi.tilebox();

//       amrex::ParallelFor(
//           bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
//             pointCFL(i, j, k, *arrayCFL, prims, lclosures, dx, dt);
//           });
//     };

//     for (int idir = 0; idir < AMREX_SPACEDIM; idir++) {
//       for (int icomp = 0; icomp < AMREX_SPACEDIM; icomp++) {
//         ParallelDescriptor::ReduceRealMax(
//             (*arrayCFL)(idir, icomp), ParallelDescriptor::IOProcessorNumber());
//       }
//     }

//     amrex::Print().SetPrecision(17) << "\n[CNS level " << level << "]\n"
//                                     << "   Rho  min, max = " << prims_min[0]
//                                     << " , " << prims_max[0] << "\n"
//                                     << "   Ux   min, max = " << prims_min[1]
//                                     << " , " << prims_max[1] << "\n"
//                                     << "   Uy   min, max = " << prims_min[2]
//                                     << " , " << prims_max[2] << "\n"
//                                     << "   Uz   min, max = " << prims_min[3]
//                                     << " , " << prims_max[3] << "\n"
//                                     << "   P    min, max = " << prims_min[5]
//                                     << " , " << prims_max[5] << "\n"
//                                     << "   T    min, max = " << prims_min[4]
//                                     << " , " << prims_max[4] << "\n";

//     amrex::Print().SetPrecision(17)
//         << "   Vax  min, max = " << (*arrayCFL)(0, 0) << " , "
//         << (*arrayCFL)(0, 1) << "\n"
//         << "   Vay  min, max = " << (*arrayCFL)(1, 0) << " , "
//         << (*arrayCFL)(1, 1) << "\n"
//         << "   Vaz  min, max = " << (*arrayCFL)(2, 0) << " , "
//         << (*arrayCFL)(2, 1) << "\n"
//         << "   CFLx         = " << (*arrayCFL)(0, 2) << "\n"
//         << "   CFLy         = " << (*arrayCFL)(1, 2) << "\n"
//         << "   CFLz         = " << (*arrayCFL)(2, 2) << "\n \n";

//     amrex::Print().SetPrecision(17) << "   Total mass   = " << tot[0] << "\n"
//                                     << "   Total x-mom  = " << tot[1] << "\n"
//                                     << "   Total y-mom  = " << tot[2] << "\n"
//                                     << "   Total z-mom  = " << tot[3] << "\n"
//                                     << "   Total energy = " << tot[4] << "\n";

// #if AMREX_USE_GPU
//     The_Arena()->free(arrayCFL);
// #else
//   delete arrayCFL;
// #endif
// #ifdef BL_LAZY
//   });
// #endif
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
  IBM::ib.ibMFa[level]->copytoRealMF(plotMF, 0, cnt);
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


Real CNS::advance(Real time, Real dt, int /*iteration*/, int /*ncycle*/) {
  BL_PROFILE("CNS::advance()");

  // Print() << "-------- advance start -------" << std::endl;
  // Print() << "time = " << time << std::endl;
  // Print() << "dt = " << dt << std::endl;
  // state[0].printTimeInterval(std::cout);
  // Print() << "---------------------------------" << std::endl;

  for (int i = 0; i < num_state_data_types; ++i) {
    state[i].allocOldData();
    state[i].swapTimeLevels(dt);
  }

  MultiFab& s1 = get_old_data(State_Type);
  MultiFab& s2 = get_new_data(State_Type);

  int ncons = d_prob_closures->NCONS;
  int nghost= d_prob_closures->NGHOST;
  MultiFab stemp(grids,dmap,ncons,nghost,MFInfo(),Factory());

  FluxRegister* fr_as_crse = nullptr;
  if (do_reflux && level < parent->finestLevel()) {
    CNS& fine_level = getLevel(level + 1);
    fr_as_crse = fine_level.flux_reg.get();
  }

  FluxRegister* fr_as_fine = nullptr;
  if (do_reflux && level > 0) {
    fr_as_fine = flux_reg.get();
  }

  if (fr_as_crse) {
    fr_as_crse->setVal(Real(0.0));
  }

  PROB::ProbLHS::advance(*this,s1,s2,stemp,fr_as_crse,fr_as_fine,time,dt,ncons,nghost);

  return dt;
}


