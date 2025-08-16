#include <AMReX_ParmParse.H>

#include <climits>
#if CNS_USE_EB
#include <AMReX_EBAmrUtil.H>
#include <AMReX_EBFArrayBox.H>
#include <AMReX_EBMultiFabUtil.H>
#endif
// #include <AMReX_PlotFileUtil.H>

#include "CNS.H"
#include "CNS_K.H"
#include "derive.H"
#include "prob.H"
#include "tagging.H"

using namespace amrex;

////////////////////////////////////////////////////////////////////////////
//    Load default parameters values                                      //
////////////////////////////////////////////////////////////////////////////

#include "default_parm.H"

amrex::Vector<amrex::Real> UniqueRand::_data;

////////////////////////////////////////////////////////////////////////////
//    Overriding AmrLevel virtual functions                               //
////////////////////////////////////////////////////////////////////////////

CNS::CNS(Amr& papa, int lev, const Geometry& level_geom, const BoxArray& bl,
         const DistributionMapping& dm, Real time)
  : AmrLevel(papa, lev, level_geom, bl, dm, time)
{
  if (do_reflux && level > 0) {
    flux_reg.define(bl, papa.boxArray(level - 1), dm,
                    papa.DistributionMap(level - 1), level_geom,
                    papa.Geom(level - 1), papa.refRatio(level - 1), level, UFA);
  }

  buildMetrics();

  Sborder.define(grids, dmap, LEN_STATE, NUM_GROW, MFInfo(), Factory());
  shock_sensor_mf.define(grids, dmap, 1, 3, MFInfo(), Factory());
  // if (!use_hybrid_scheme) {
  //   shock_sensor_mf.setVal(1.0); // default to shock-capturing scheme
  // }
  // ifine_mask.define(grids, dmap, 1, 0, MFInfo());
  // fillFineMask();

  // Initialize reaction
  get_new_data(Reactions_Type).setVal(0.0);
  if (do_react) {
    if (chem_integrator == "ReactorNull") {
      Print() << "WARNING: turning on reactions while using ReactorNull. "
                 "Make sure this is intended."
              << std::endl;
    }

    reactor = pele::physics::reactions::ReactorBase::create(chem_integrator);
    reactor->init(1, 1); // init reactor to constant volume type
  }

#if (NUM_FIELD > 0)
  // Initialise random number
  Vector<IntVect> ref_ratio;
  for (int lv = 0; lv < parent->maxLevel(); ++lv) {
    ref_ratio.push_back(parent->refRatio(lv));
  }
  WienerProcess.init(AMREX_SPACEDIM, level, ref_ratio);
#endif

  // Initialise LES model
  if (do_les || do_pasr) {
    ParmParse pp("cns");
    std::string les_model_name;
    pp.get("les_model", les_model_name);
    les_model = LESModel::create(les_model_name);
  }
}

CNS::~CNS()
{
  if (do_react) { reactor->close(); }
}

void CNS::init(AmrLevel& old)
{
  BL_PROFILE("CNS::init(old)");

  auto& oldlev = dynamic_cast<CNS&>(old);

  Real dt_new = parent->dtLevel(level);
  Real cur_time = oldlev.state[State_Type].curTime();
  Real prev_time = oldlev.state[State_Type].prevTime();
  Real dt_old = cur_time - prev_time;
  setTimeLevel(cur_time, dt_old, dt_new);

  MultiFab& S_new = get_new_data(State_Type);
  FillPatch(old, S_new, 0, cur_time, State_Type, 0, LEN_STATE);

  amrex::MultiFab& React_new = get_new_data(Reactions_Type);
  if (do_react) {
    FillPatch(old, React_new, 0, cur_time, Reactions_Type, 0, React_new.nComp());
  } else {
    React_new.setVal(0.0);
  }

  MultiFab& C_new = get_new_data(Cost_Type);
  FillPatch(old, C_new, 0, cur_time, Cost_Type, 0, 1);
}

void CNS::init()
{
  // This version inits the data on a new level that did not
  // exist before regridding.
  BL_PROFILE("CNS::init()");

  Real dt = parent->dtLevel(level);
  Real cur_time = getLevel(level - 1).state[State_Type].curTime();
  Real prev_time = getLevel(level - 1).state[State_Type].prevTime();
  Real dt_old =
    (cur_time - prev_time) / static_cast<Real>(parent->MaxRefRatio(level - 1));
  setTimeLevel(cur_time, dt_old, dt);

  MultiFab& S_new = get_new_data(State_Type);
  FillCoarsePatch(S_new, 0, cur_time, State_Type, 0, LEN_STATE);

  MultiFab& C_new = get_new_data(Cost_Type);
  FillCoarsePatch(C_new, 0, cur_time, Cost_Type, 0, 1);
}

void CNS::initData()
{
  BL_PROFILE("CNS::initData()");

  if (verbose > 0) {
    Print() << "Initializing the data at level " << level << std::endl;
  }

  MultiFab& S_new = get_new_data(State_Type);
  S_new.setVal(0.0);

  // <T> const* = pointer to constant <T>; const <T>* == <T> const*
  ProbParm const* lprobparm = d_prob_parm;
  const auto geomdata = geom.data();
#ifdef USE_PMFDATA
  pele::physics::PMF::PmfData::DataContainer const* lpmfdata =
    pmf_data.getDeviceData();
#endif

  auto const& sarrs = S_new.arrays();
  amrex::ParallelFor(S_new,
                     [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k) noexcept {
                       prob_initdata(i, j, k, sarrs[box_no], geomdata, *lprobparm
#ifdef USE_PMFDATA
                                     ,
                                     lpmfdata
#endif
                       );

                       // Verify that the sum of (rho Y)_i = rho at every cell
                       cns_check_species_sum_to_one(i, j, k, sarrs[box_no]);
                     });
  enforce_consistent_state(S_new);

  get_new_data(Reactions_Type).setVal(0.0);
  get_new_data(Cost_Type).setVal(1.0);
}

void CNS::computeInitialDt(int finest_level, int /*sub_cycle*/, Vector<int>& n_cycle,
                           const Vector<IntVect>& /*ref_ratio*/,
                           Vector<Real>& dt_level, Real stop_time)
{
  // Grids have been constructed, compute dt for all levels
  if (level > 0) return;

  // Find the minimum over all levels
  Real dt_0 = std::numeric_limits<Real>::max();
  int n_factor = 1;
  for (int i = 0; i <= finest_level; i++) {
    dt_level[i] = getLevel(i).estTimeStep();
    n_factor *= n_cycle[i];
    dt_0 = std::min<amrex::Real>(dt_0, n_factor * dt_level[i]);
  }

  // Limit dt's by the value of stop_time
  const Real eps = 0.001 * dt_0;
  Real cur_time = state[State_Type].curTime();
  if (stop_time >= 0.0) {
    if ((cur_time + dt_0) > (stop_time - eps)) dt_0 = stop_time - cur_time;
  }

  n_factor = 1;
  for (int i = 0; i <= finest_level; i++) {
    n_factor *= n_cycle[i];
    dt_level[i] = dt_0 / n_factor;
  }
}

void CNS::computeNewDt(int finest_level, int /*sub_cycle*/, Vector<int>& n_cycle,
                       const Vector<IntVect>& /*ref_ratio*/, Vector<Real>& dt_min,
                       Vector<Real>& dt_level, Real stop_time, int post_regrid_flag)
{
  // We are at the end of a coarse grid timecycle.
  // Compute the timesteps for the next iteration.
  if (level > 0) return;

  for (int i = 0; i <= finest_level; i++) { dt_min[i] = getLevel(i).estTimeStep(); }

  if (post_regrid_flag == 1) {
    // Limit dt's by pre-regrid dt
    for (int i = 0; i <= finest_level; i++) {
      dt_min[i] = std::min(dt_min[i], dt_level[i]);
    }
  } else {
    // Limit dt's by change_max * old dt (TODO: Do we really need this?)
    for (int i = 0; i <= finest_level; i++) {
      if ((verbose > 0) && ParallelDescriptor::IOProcessor()) {
        Print() << " >> Compute new dt: limiting dt at level " << i << '\n';
        Print() << "    new dt computed: " << dt_min[i] << '\n';
        if (dt_min[i] > dt_max_change * dt_level[i]) {
          Print() << "    but limited to: " << dt_max_change * dt_level[i] << " = "
                  << dt_max_change << " * " << dt_level[i] << '\n';
        }
      }

      dt_min[i] = std::min(dt_min[i], dt_max_change * dt_level[i]);
    }
  }

  // Find the minimum over all levels
  Real dt_0 = std::numeric_limits<Real>::max();
  int n_factor = 1;
  for (int i = 0; i <= finest_level; i++) {
    n_factor *= n_cycle[i];
    dt_0 = std::min(dt_0, n_factor * dt_min[i]);
  }

  // Limit dt's by the value of stop_time.
  const Real eps = 0.001 * dt_0;
  Real cur_time = state[State_Type].curTime();
  if (stop_time >= 0.0) {
    if ((cur_time + dt_0) > (stop_time - eps)) { dt_0 = stop_time - cur_time; }
  }

  n_factor = 1;
  for (int i = 0; i <= finest_level; i++) {
    n_factor *= n_cycle[i];
    dt_level[i] = dt_0 / n_factor;
  }
}

void CNS::post_regrid(int /*lbase*/, int /*new_finest*/)
{
  enforce_consistent_state();

  if (do_react && use_typical_vals_chem) { set_typical_values_chem(); }

  // fillFineMask();
}

void CNS::post_timestep(int iteration)
{
  BL_PROFILE("CNS::post_timestep()");

  if (do_reflux && level < parent->finestLevel()) {
    if (verbose > 0) { amrex::Print() << " >> Reflux() at level " << level << '\n'; }

    CNS& fine_level = getLevel(level + 1);
    MultiFab& S_crse = get_new_data(State_Type);
    const int ncomp = UFA;
#if CNS_USE_EB
    MultiFab& S_fine = fine_level.get_new_data(State_Type);
    fine_level.flux_reg.Reflux(S_crse, *volfrac, S_fine, *fine_level.volfrac, 0, 0,
                               ncomp);
#else
    fine_level.flux_reg.Reflux(S_crse, 0, 0, ncomp);
#endif
    enforce_consistent_state(S_crse);
  }

  if (level < parent->finestLevel()) {
    avgDown();
    getLevel(level + 1).resetFillPatcher();
  }

#ifdef USE_FULL_PROB_POST_TIMESTEP
  full_prob_post_timestep(iteration);
#else
  MultiFab& S = get_new_data(State_Type);
  MultiFab& IR = get_new_data(Reactions_Type);
  Real curtime = state[State_Type].curTime();
  Real dtlev = parent->dtLevel(level);
  ProbParm const* lprobparm = d_prob_parm;
  const auto geomdata = geom.data();
  auto const& sarrs = S.arrays();
  auto const& irarrs = IR.const_arrays();
  amrex::ParallelFor(S,
                     [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k) noexcept {
                       // Do sth here, e.g. calculate time statistics
                       prob_post_timestep(i, j, k, curtime, dtlev, sarrs[box_no],
                                          irarrs[box_no], geomdata, *lprobparm);
                     });
#endif

  if (do_react && use_typical_vals_chem &&
      parent->levelSteps(0) % reset_typical_vals_int == 0) {
    set_typical_values_chem();
  }
}

void CNS::postCoarseTimeStep(Real time)
{
  BL_PROFILE("CNS::postCoarseTimeStep()");

  AmrLevel::postCoarseTimeStep(time);

#ifdef USE_PROB_POST_COARSETIMESTEP
  prob_post_coarsetimestep(time);
#endif

  printTotalandCheckNan(); // must do because this checks for nan as well

  checkRuntimeMessages(); // re-read inputs or soft exit simulation
}

void CNS::post_init(Real /*stop_time*/)
{
  // Initialise reaction
  if (do_react) {
    if (use_typical_vals_chem) { set_typical_values_chem(); }

    react_state(parent->cumTime(), parent->dtLevel(level), true);
  }

  // {
  //   amrex::MultiFab ifine_mask_plot(grids, dmap, 1, 0);
  //   auto const& ifm_arrs = ifine_mask.const_arrays();
  //   auto const& ifmp_arrs = ifine_mask_plot.arrays();
  //   amrex::ParallelFor(ifine_mask_plot,
  //     [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k) noexcept {
  //       ifmp_arrs[box_no](i, j, k) = ifm_arrs[box_no](i, j, k);
  //     });
  //   amrex::WriteSingleLevelPlotfile("plt_lev" + std::to_string(level),
  //   ifine_mask_plot, {"mask"}, geom, 0.0, 0);
  // }

  if (level > 0) return;

  // Average data down from finer levels
  // so that conserved data is consistent between levels
  for (int k = parent->finestLevel() - 1; k >= 0; --k) { getLevel(k).avgDown(); }

  printTotalandCheckNan();
}

void CNS::post_restart()
{
  // Copy problem parameter structs to device
  // amrex::Gpu::copy(amrex::Gpu::hostToDevice, h_prob_parm,
  //                  h_prob_parm+1, d_prob_parm);
  // Done in prob.H

#ifdef AMREX_USE_GPU
  // Cannot use Gpu::copy because ProbParm is not trivailly copyable.
  Gpu::htod_memcpy_async(CNS::d_prob_parm, CNS::h_prob_parm, sizeof(ProbParm));
#else
  std::memcpy(CNS::d_prob_parm, CNS::h_prob_parm, sizeof(ProbParm));
#endif

  Sborder.define(grids, dmap, LEN_STATE, NUM_GROW, MFInfo(), Factory());
  shock_sensor_mf.define(grids, dmap, 1, 3, MFInfo(), Factory());
  // if (!use_hybrid_scheme) {
  //   shock_sensor_mf.setVal(1.0); // default to shock-capturing scheme
  // }
  // ifine_mask.define(grids, dmap, 1, 0, MFInfo());
  // fillFineMask();

  MultiFab& S_new = get_new_data(State_Type);

#if (NUM_FIELD > 0)
  // Initialise random number
  Vector<IntVect> ref_ratio;
  for (int lv = 0; lv < parent->maxLevel(); ++lv) {
    ref_ratio.push_back(parent->refRatio(lv));
  }
  WienerProcess.init(AMREX_SPACEDIM, level, ref_ratio);

  // Populate fields (when restarting from a different number of fields)
  if (do_restart_fields) {
    Print() << " >> Resetting stochastic fields state data ..." << std::endl;

    // Move aux variables
    MultiFab::Copy(S_new, S_new, NVAR, UFA, NUM_AUX, 0);

    // Copy mean to fields
    for (int nf = 1; nf < NUM_FIELD + 1; ++nf) {
      MultiFab::Copy(S_new, S_new, 0, nf * NVAR, NVAR, 0);
    }

    if (do_react) {
      Print() << " >> Resetting stochastic fields reaction data ..." << std::endl;

      MultiFab& I_R = get_new_data(Reactions_Type);

      // Copy mean to fields
      for (int nf = 1; nf < NUM_FIELD + 1; ++nf) {
        MultiFab::Copy(I_R, I_R, 0, nf * NREACT, NREACT, 0);
      }
    }
  }
#endif

  // Initialize reactor
  if (do_react) {
    if (chem_integrator == "ReactorNull") {
      amrex::Print() << "WARNING: turning on reactions while using ReactorNull. "
                        "Make sure this is intended.\n";
    }

    reactor = pele::physics::reactions::ReactorBase::create(chem_integrator);
    reactor->init(1, 1);

    if (use_typical_vals_chem) { set_typical_values_chem(); }

    react_state(parent->cumTime(), parent->dtLevel(level), true);
  }

  // Initialise LES model
  if (do_les || do_pasr) {
    ParmParse pp("cns");
    std::string les_model_name;
    pp.get("les_model", les_model_name);
    les_model = LESModel::create(les_model_name);
  }

  ProbParm const* lprobparm = d_prob_parm;
  const auto geomdata = geom.data();
  auto const& sarrs = S_new.arrays();
  amrex::ParallelFor(
    S_new, [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k) noexcept {
      // Modify restarted data and/or add turbulence
      prob_post_restart(i, j, k, sarrs[box_no], geomdata, *lprobparm);
    });
  enforce_consistent_state(S_new);
}

void CNS::errorEst(TagBoxArray& tags, int clearval, int tagval, Real time,
                   int /*n_error_buf = 0*/, int /*ngrow*/)
{
  BL_PROFILE("CNS::errorEst()");

#if CNS_USE_EB
  if (refine_cutcells && (level < refine_cutcells_max_lev)) {
    const MultiFab& S_new = get_new_data(State_Type);
    TagCutCells(tags, S_new);
  }
#endif

  for (int j = 0; j < errtags.size(); ++j) {
    std::unique_ptr<MultiFab> mf;
    if (errtags[j].Field() != std::string()) {
      mf = derive(errtags[j].Field(), time, errtags[j].NGrow());
    }
    errtags[j](tags, mf.get(), clearval, tagval, time, level, geom);
  }

  // Tag user specified box(es)
  if (!refine_boxes.empty()) {
    const int n_refine_boxes = refine_boxes.size();
    const auto problo = geom.ProbLoArray();
    const auto dx = geom.CellSizeArray();
    auto boxes = dp_refine_boxes;

    auto const& tagma = tags.arrays();
    ParallelFor(
      tags, [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k) noexcept {
        RealVect pos{AMREX_D_DECL((i + 0.5) * dx[0] + problo[0],
                                  (j + 0.5) * dx[1] + problo[1],
                                  (k + 0.5) * dx[2] + problo[2])};
        for (int irb = 0; irb < n_refine_boxes; ++irb) {
          if (boxes[irb].contains(pos) && (level < refine_boxes_max_lev[irb])) {
            tagma[box_no](i, j, k) = tagval;
          }
        }
      });
  }

  const MultiFab& S_new = get_new_data(State_Type);
  MultiFab Sg1(grids, dmap, NVAR, 1, MFInfo(), Factory());
  const Real cur_time = state[State_Type].curTime();
  FillPatch(*this, Sg1, Sg1.nGrow(), cur_time, State_Type, URHO, NVAR, 0);

  auto const& tagma = tags.arrays();
  const auto dxinv = geom.InvCellSizeArray();

  if (level < refine_dengrad_max_lev) {
    const Real threshold = refine_dengrad[level];

    MultiFab rho(Sg1, amrex::make_alias, URHO, 1);
    auto const& rhoma = rho.const_arrays();

    ParallelFor(rho, IntVect(-1), [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k) noexcept {
      cns_tag_grad(i, j, k, tagma[box_no], rhoma[box_no], dxinv, threshold, tagval);
    });
  }

  if (level < refine_velgrad_max_lev) {
    const Real threshold = refine_velgrad[level];

    for (MFIter mfi(tags, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
      const Box& bx = mfi.tilebox();

      FArrayBox velgrad(bx, 1);
      const int* bc = 0; // dummy
      cns_dervelgrad(bx, velgrad, 0, 1, Sg1[mfi], geom, cur_time, bc, level);

      auto tagarr = tags.array(mfi);
      auto mvarr = velgrad.const_array();
      ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        cns_tag_error(i, j, k, tagarr, mvarr, threshold, tagval);
      });
    }
  }

  if (level < refine_presgrad_max_lev) {
    const Real threshold = refine_presgrad[level];

    for (MFIter mfi(tags, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
      const Box& bx = mfi.tilebox();

      const Box& bxg1 = amrex::grow(bx, 1);
      FArrayBox pres(bxg1, 1);
      const int* bc = 0; // dummy
      cns_derpres(bxg1, pres, 0, 1, Sg1[mfi], geom, cur_time, bc, level);

      auto tagarr = tags.array(mfi);
      auto parr = pres.const_array();
      ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        cns_tag_grad(i, j, k, tagarr, parr, dxinv, threshold, tagval);
      });
    }
  }

  if (level < refine_magvort_max_lev) {
    const Real threshold = refine_magvort[level];

    for (MFIter mfi(tags, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
      const Box& bx = mfi.tilebox();

      FArrayBox magvort(bx, 1);
      const int* bc = 0; // dummy
      cns_dermagvort(bx, magvort, 0, 1, Sg1[mfi], geom, cur_time, bc, level);

      auto tagarr = tags.array(mfi);
      auto mvarr = magvort.const_array();
      ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        cns_tag_error(i, j, k, tagarr, mvarr, threshold, tagval);
      });
    }
  }

#if (NUM_FIELD > 0) // cannot calculate tke if no SF
  if (level < refine_tke_max_lev) {
    const Real threshold = refine_tke[level];

    for (MFIter mfi(tags, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
      const Box& bx = mfi.tilebox();

      FArrayBox tke(bx, 1);
      const int* bc = 0; // dummy
      cns_dertke(bx, tke, 0, 1, S_new[mfi], geom, cur_time, bc, level);

      auto tagarr = tags.array(mfi);
      auto karr = tke.const_array();
      ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        cns_tag_error(i, j, k, tagarr, karr, threshold, tagval);
      });
    }
  }
#endif

  // Problem specific taggings
  ProbParm const* lprobparm = d_prob_parm;
  const auto geomdata = geom.data();
  auto tagarr = tags.arrays();
  auto const& sarrs = S_new.arrays();
  amrex::ParallelFor(tags,
                     [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k) noexcept {
                       prob_tag_error(i, j, k, tagarr[box_no], sarrs[box_no], level,
                                      tagval, time, geomdata, *lprobparm);
                     });

  Gpu::streamSynchronize();
}

////////////////////////////////////////////////////////////////////////////
//    Helper functions that are in the CNS class but not in AMRLevel      //
////////////////////////////////////////////////////////////////////////////

/**
 * \brief Print total and check for nan
 */
void CNS::printTotalandCheckNan() const
{
  BL_PROFILE("CNS::printTotalandCheckNan()");

  const MultiFab& S_new = get_new_data(State_Type);
  MultiFab mf(grids, dmap, 1, 0);
  std::array<Real, 6> tot;
  for (int comp = 0; comp < 6; ++comp) {
    MultiFab::Copy(mf, S_new, comp, 0, 1, 0);
#if CNS_USE_EB
    MultiFab::Multiply(mf, *volfrac, 0, 0, 1, 0);
#endif
    tot[comp] = mf.sum(0, true) * geom.ProbSize();
  }
#ifdef BL_LAZY
  Lazy::QueueReduction([=]() mutable {
#endif
    ParallelDescriptor::ReduceRealSum(tot.data(), 5,
                                      ParallelDescriptor::IOProcessorNumber());
    if (verbose > 0) {
      amrex::Print().SetPrecision(17)
        << "\n[CNS] Total mass       = " << tot[0]
        << "\n" AMREX_D_TERM(<< "      Total x-momentum = " << tot[1] << "\n",
                             << "      Total y-momentum = " << tot[2] << "\n",
                             << "      Total z-momentum = " << tot[3] << "\n")
        << "      Total energy     = " << tot[4] << "\n";
    }
#ifdef BL_LAZY
  });
#endif

  // Nan detector for soft exit
  if (isnan(tot[0]) ||
      AMREX_D_TERM(isnan(tot[1]), || isnan(tot[2]), || isnan(tot[3])) ||
      isnan(tot[4])) {
    signalStopJob = true;
  }
  amrex::ParallelDescriptor::ReduceBoolOr(signalStopJob);
}

/**
 * \brief Allow re-reading inputs file or soft exit of simulation at runtime
 */
void CNS::checkRuntimeMessages() 
{
  if (parent->levelSteps(0) % check_message_int == 0) {    
    int action_flag = 0;
    if (ParallelDescriptor::IOProcessor()) {
      std::string action_file = "reread_inputs";
      FILE* fp = fopen(action_file.c_str(), "r");
      if (fp != nullptr) {
        remove(action_file.c_str());
        action_flag = 1;
        fclose(fp);
      }

      action_file = "stop_simulation";
      fp = fopen(action_file.c_str(), "r");
      if (fp != nullptr) {
        remove(action_file.c_str());
        action_flag = 2;
        fclose(fp);
      }
    }

    amrex::ParallelDescriptor::Bcast(
      &action_flag, 1, ParallelDescriptor::IOProcessorNumber());

    // if (action_flag == 1) {
    //   amrex::Print() << " >> Re-reading inputs ...\n";
    //   for (int lev = 0; lev <= parent->finestLevel(); ++lev) {
    //     getLevel(lev).read_params(); // Cannot do, need to reinit paramparse
    //   }
    // } else 
    if (action_flag == 2) {
      amrex::Print() << " >> Stopping simulation ...\n";
      signalStopJob = true;
    }
  }
}

// /**
//  * \brief Write time statistics to a file
//  */
// void
// CNS::writeTimeStat ()
// {
//   BL_PROFILE("CNS::writeTimeStat()");

//   const Real dtlev = parent->dtLevel(level);

//   MultiFab& S_data = get_new_data(State_Type);
//   MultiFab& I_R = get_new_data(Reactions_Type);

// #ifdef AMREX_USE_OMP
// #pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
// #endif
// {
//   for (MFIter mfi(S_data, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
//     const Box& tilebox = mfi.tilebox();
//     auto Sarr = S_data.array(mfi);
//     const auto IRarr = I_R.array(mfi);

//     amrex::ParallelFor(tilebox,
//     [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
//     });
//   }
// }
// }

/**
 * \brief Return a flag signalling is it okay to continue the simulation.
 *        test == 1: ok, test == 0: stop.
 */
int CNS::okToContinue()
{
  if (level > 0) { return 1; }

  int test = 1;
  if (signalStopJob) {
    test = 0;
    amrex::Print() << "Signalling a stop of the run due to signalStopJob = true."
                   << std::endl;
  } else if (parent->dtLevel(0) < dt_cutoff) {
    test = 0;
    amrex::Print() << "Signalling a stop of the run because dt < dt_cutoff."
                   << std::endl;
  }

  return test;
}

void CNS::avgDown()
{
  BL_PROFILE("CNS::avgDown()");

  if (level == parent->finestLevel()) return;

  auto& fine_lev = getLevel(level + 1);
  MultiFab& S_crse = get_new_data(State_Type);
  const MultiFab& S_fine = fine_lev.get_new_data(State_Type);

#if CNS_USE_EB
  MultiFab volume(S_fine.boxArray(), S_fine.DistributionMap(), 1, 0);
  geom.GetVolume(volume);
  amrex::EB_average_down(S_fine, S_crse, volume, fine_lev.volFrac(), 0,
                         S_fine.nComp(), fine_ratio);
#else
  const amrex::Geometry& fgeom = getLevel(level + 1).geom;
  const amrex::Geometry& cgeom = geom;
  amrex::average_down(S_fine, S_crse, fgeom, cgeom, 0, S_fine.nComp(), fine_ratio);
#endif

  if (do_react) {
    MultiFab& R_crse = get_new_data(Reactions_Type);
    const MultiFab& R_fine = fine_lev.get_new_data(Reactions_Type);
#if CNS_USE_EB
    amrex::EB_average_down(R_fine, R_crse, volume, fine_lev.volFrac(), 0,
                           R_fine.nComp(), fine_ratio);
#else
    amrex::average_down(R_fine, R_crse, fgeom, cgeom, 0, R_fine.nComp(), fine_ratio);
#endif
  }
}

void CNS::buildMetrics()
{
  BL_PROFILE("CNS::buildMetrics()");

#if CNS_USE_EB
  static_assert(AMREX_SPACEDIM > 1, "EB only supports 2D and 3D");

  ParmParse ppeb2("eb2");
  std::string geom_type = "all_regular";
  ppeb2.query("geom_type", geom_type);

  if (geom_type != "all_regular") {
    // make sure dx == dy == dz if use EB
    const Real* dx = geom.CellSize();
    if (AMREX_D_TERM(, std::abs(dx[0] - dx[1]) > 1.e-12 * dx[0],
                     || std::abs(dx[0] - dx[2]) > 1.e-12 * dx[0])) {
      amrex::Abort("EB must have dx == dy == dz (for cut surface fluxes)\n");
    }
  }

  const auto& ebfactory = dynamic_cast<EBFArrayBoxFactory const&>(Factory());
  volfrac = &(ebfactory.getVolFrac());
  bndrycent = &(ebfactory.getBndryCent());
  areafrac = ebfactory.getAreaFrac();
  facecent = ebfactory.getFaceCent();

  // Level mask for redistribution
  level_mask.clear();
  level_mask.define(grids, dmap, 1, 3);
  level_mask.BuildMask(
    geom.Domain(), geom.periodicity(), CNSConstants::level_mask_covered,
    CNSConstants::level_mask_notcovered, CNSConstants::level_mask_physbnd,
    CNSConstants::level_mask_interior);
#endif
}

Real CNS::estTimeStep()
{
  BL_PROFILE("CNS::estTimeStep()");

  if (fixed_dt > 0) { return fixed_dt; }

  const auto dx = geom.CellSizeArray();
  const MultiFab& S = get_new_data(State_Type);

#if CNS_USE_EB
  auto const& fact = dynamic_cast<EBFArrayBoxFactory const&>(S.Factory());
  auto const& flags = fact.getMultiEBCellFlagFab();
#endif

  Real estdt = std::numeric_limits<Real>::max();

  // Reduce min operation
  ReduceOps<ReduceOpMin> reduce_op;
  ReduceData<Real> reduce_data(reduce_op);
  using ReduceTuple = typename decltype(reduce_data)::Type;

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
  for (MFIter mfi(S, false); mfi.isValid(); ++mfi) {
    const Box& bx = mfi.tilebox();
    auto const& s_arr = S.array(mfi);
#if CNS_USE_EB
    const auto& flag = flags[mfi];
    if (flag.getType(bx) != FabType::covered)
#endif
    {
      reduce_op.eval(bx, reduce_data,
                     [=] AMREX_GPU_DEVICE(int i, int j, int k) -> ReduceTuple {
                       return cns_estdt_hydro(i, j, k, s_arr, dx);
                     });

      if (do_visc) {
        auto const* ltransparm = trans_parms.device_trans_parm();
        reduce_op.eval(bx, reduce_data,
                       [=] AMREX_GPU_DEVICE(int i, int j, int k) -> ReduceTuple {
                         return cns_estdt_visc(i, j, k, s_arr, dx, ltransparm);
                       });
      }
    } // if not covered
  } // mfi

  ReduceTuple host_tuple = reduce_data.value();
  estdt = amrex::min(estdt, amrex::get<0>(host_tuple));

  estdt *= cfl;
  ParallelDescriptor::ReduceRealMin(estdt);
  return estdt;
}

// void CNS::fillFineMask()
// {
//   BL_PROFILE("CNS::fillFineMask()");

//   if (level == parent->finestLevel()) {
//     ifine_mask.setVal(1);
//     return;
//   }

//   iMultiFab tmp_mask(grids, dmap, 1, 2, MFInfo());
//   tmp_mask = makeFineMask(tmp_mask, parent->boxArray(level + 1), fine_ratio,
//   geom.periodicity(), 1, 0);

//   // Grow by 2 cells
//   auto const& tmask_arrs = tmp_mask.const_arrays();
//   auto const& fmask_arrs = ifine_mask.arrays();
//   amrex::ParallelFor(ifine_mask,
//     [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k) noexcept {
//       fmask_arrs[box_no](i, j, k) = 0;
//       for (int ii = -2; ii <= 2; ++ii) {
//         for (int jj = -2; jj <= 2; ++jj) {
//           for (int kk = -2; kk <= 2; ++kk) {
//             if (tmask_arrs[box_no](i + ii, j + jj, k + kk) == 1) {
//               fmask_arrs[box_no](i, j, k) = 1;
//             }
//           }
//         }
//       }
//     });
// }