#include <climits>

#include <AMReX_ParmParse.H>
// #include <AMReX_ErrorList.H>
#if CNS_USE_EB
#include <AMReX_EBMultiFabUtil.H>
#include <AMReX_EBFArrayBox.H>
#include <AMReX_EBAmrUtil.H>
#endif

#include "CNS.H"
#include "CNS_K.H"
#include "tagging.H"
#include "derive.H"
#include "prob.H"

using namespace amrex;

////////////////////////////////////////////////////////////////////////////
//    Load default parameters values                                      //
////////////////////////////////////////////////////////////////////////////

#include "default_parm.H"

amrex::Vector<amrex::Real> UniqueRand::_data;

////////////////////////////////////////////////////////////////////////////
//    Overriding AmrLevel virtual functions                               //
////////////////////////////////////////////////////////////////////////////

CNS::CNS (Amr&            papa,
          int             lev,
          const Geometry& level_geom,
          const BoxArray& bl,
          const DistributionMapping& dm,
          Real            time) : AmrLevel(papa,lev,level_geom,bl,dm,time)
{
  if (do_reflux && level > 0) {
    flux_reg.define(bl, papa.boxArray(level-1),
                    dm, papa.DistributionMap(level-1),
                    level_geom, papa.Geom(level-1),
                    papa.refRatio(level-1), level, LEN_STATE);
  }

  buildMetrics();

  // Initialize reaction
  get_new_data(Reactions_Type).setVal(0.0);
  if (do_react) {
    if (chem_integrator == "ReactorNull") {
      Print() << "WARNING: turning on reactions while using ReactorNull. "
                 "Make sure this is intended." << std::endl;
    }

    reactor = pele::physics::reactions::ReactorBase::create(chem_integrator);
    reactor->init(1, 1); //init reactor to constant volume type
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
    les_model = LESModel::create(les_model_name);
  }

#if CNS_USE_EB
  // make sure dx = dy = dz
  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx = geom.CellSizeArray();
  const amrex::Real small = 1.e-13;
  if (amrex::max<amrex::Real>(AMREX_D_DECL(
      static_cast<amrex::Real>(0.0),
      static_cast<amrex::Real>(std::abs(dx[0] - dx[1])),
      static_cast<amrex::Real>(std::abs(dx[0] - dx[2])))) > small * dx[0]) {
    amrex::Abort("EB formulation assumes dx == dy == dz");
  }
#endif
}

CNS::~CNS () 
{
  if (do_react) {
    reactor->close();
  }
}

void
CNS::init (AmrLevel& old)
{
  BL_PROFILE("CNS::init(old)");

  auto& oldlev = dynamic_cast<CNS&>(old);

  Real dt_new    = parent->dtLevel(level);
  Real cur_time  = oldlev.state[State_Type].curTime();
  Real prev_time = oldlev.state[State_Type].prevTime();
  Real dt_old    = cur_time - prev_time;
  setTimeLevel(cur_time,dt_old,dt_new);

  MultiFab& S_new = get_new_data(State_Type);
  FillPatch(old, S_new, 0, cur_time, State_Type, 0, LEN_STATE);

  amrex::MultiFab& React_new = get_new_data(Reactions_Type);
  if (do_react) {
    FillPatch(old, React_new, 0, cur_time, Reactions_Type, 0, React_new.nComp());
  } else {
    React_new.setVal(0);
  }

  MultiFab& C_new = get_new_data(Cost_Type);
  FillPatch(old, C_new, 0, cur_time, Cost_Type, 0, 1);
}

void
CNS::init ()
{
  // This version inits the data on a new level that did not
  // exist before regridding.
  BL_PROFILE("CNS::init()");

  Real dt        = parent->dtLevel(level);
  Real cur_time  = getLevel(level-1).state[State_Type].curTime();
  Real prev_time = getLevel(level-1).state[State_Type].prevTime();
  Real dt_old = (cur_time - prev_time)/static_cast<Real>(parent->MaxRefRatio(level-1));
  setTimeLevel(cur_time,dt_old,dt);

  MultiFab& S_new = get_new_data(State_Type);
  FillCoarsePatch(S_new, 0, cur_time, State_Type, 0, LEN_STATE);

  MultiFab& C_new = get_new_data(Cost_Type);
  FillCoarsePatch(C_new, 0, cur_time, Cost_Type, 0, 1);
}

void
CNS::initData ()
{
  BL_PROFILE("CNS::initData()");

  if (verbose > 0) {
    Print() << "Initializing the data at level " << level << std::endl;
  }

  MultiFab& S_new = get_new_data(State_Type);
  S_new.setVal(0.0);

  Parm const* lparm = d_parm;
  ProbParm const* lprobparm = d_prob_parm; // <T> const* = pointer to constant <T>; const <T>* == <T> const*
  const auto geomdata = geom.data();

  auto const& sarrs = S_new.arrays();
  amrex::ParallelFor(S_new,
  [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept {
    prob_initdata(i, j, k, sarrs[box_no], geomdata, *lparm, *lprobparm);
    cns_check_species_sum_to_one(i, j, k, sarrs[box_no]); // Verify that the sum of (rho Y)_i = rho at every cell
  });
  amrex::Gpu::synchronize();

  get_new_data(Reactions_Type).setVal(0.0);

  get_new_data(Cost_Type).setVal(1.0);
}

void
CNS::computeInitialDt (int                    finest_level,
                       int                    /*sub_cycle*/,
                       Vector<int>&           n_cycle,
                       const Vector<IntVect>& /*ref_ratio*/,
                       Vector<Real>&          dt_level,
                       Real                   stop_time)
{
  // Grids have been constructed, compute dt for all levels
  if (level > 0) return;

  // Find the minimum over all levels
  Real dt_0 = std::numeric_limits<Real>::max();
  int n_factor = 1;
  for (int i = 0; i <= finest_level; i++) {
    dt_level[i] = getLevel(i).estTimeStep();
    n_factor *= n_cycle[i];
    dt_0 = std::min<amrex::Real>(dt_0, n_factor*dt_level[i]);
  }

  // Limit dt's by the value of stop_time
  const Real eps = 0.001*dt_0;
  Real cur_time  = state[State_Type].curTime();
  if (stop_time >= 0.0) {
    if ((cur_time + dt_0) > (stop_time - eps))
      dt_0 = stop_time - cur_time;
  }

  n_factor = 1;
  for (int i = 0; i <= finest_level; i++) {
    n_factor *= n_cycle[i];
    dt_level[i] = dt_0/n_factor;
  }
}

void
CNS::computeNewDt (int                    finest_level,
                   int                    /*sub_cycle*/,
                   Vector<int>&           n_cycle,
                   const Vector<IntVect>& /*ref_ratio*/,
                   Vector<Real>&          dt_min,
                   Vector<Real>&          dt_level,
                   Real                   stop_time,
                   int                    post_regrid_flag)
{
  // We are at the end of a coarse grid timecycle.
  // Compute the timesteps for the next iteration.
  if (level > 0) return;

  for (int i = 0; i <= finest_level; i++) {
    dt_min[i] = getLevel(i).estTimeStep();
  }

  if (post_regrid_flag == 1) {
    // Limit dt's by pre-regrid dt
    for (int i = 0; i <= finest_level; i++) {
      dt_min[i] = std::min(dt_min[i], dt_level[i]);
    }
  } else {
    // Limit dt's by change_max * old dt (TODO: Do we really need this?)
    static Real change_max = 1.1;
    for (int i = 0; i <= finest_level; i++) {
      if ((verbose > 0) && ParallelDescriptor::IOProcessor()) {                
        Print() << " >> Compute new dt: limiting dt at level " << i << '\n';
        Print() << "    new dt computed: " << dt_min[i] << '\n';
        if (dt_min[i] > change_max * dt_level[i]) {
          Print() << "    but limited to: " << change_max * dt_level[i] 
                  << " = " << change_max << " * " << dt_level[i] << '\n';
        }
      }

      dt_min[i] = std::min(dt_min[i], change_max*dt_level[i]);
    }
  }

  // Find the minimum over all levels
  Real dt_0 = std::numeric_limits<Real>::max();
  int n_factor = 1;
  for (int i = 0; i <= finest_level; i++) {
    n_factor *= n_cycle[i];
    dt_0 = std::min(dt_0, n_factor*dt_min[i]);
  }

  // Limit dt's by the value of stop_time.
  const Real eps = 0.001*dt_0;
  Real cur_time  = state[State_Type].curTime();
  if (stop_time >= 0.0) {
    if ((cur_time + dt_0) > (stop_time - eps)) {
      dt_0 = stop_time - cur_time;
    }
  }

  n_factor = 1;
  for (int i = 0; i <= finest_level; i++) {
    n_factor *= n_cycle[i];
    dt_level[i] = dt_0/n_factor;
  }
}

void
CNS::post_regrid (int /*lbase*/, int /*new_finest*/)
{
  if (do_react && use_typical_vals_chem) {
    set_typical_values_chem();        
  }
  
  enforce_consistent_state();
}

void
CNS::post_timestep (int /*iteration*/)
{
  BL_PROFILE("CNS::post_timestep()");

  if (do_reflux && level < parent->finestLevel()) {
    CNS& fine_level = getLevel(level+1);
    MultiFab& S_crse = get_new_data(State_Type);
    MultiFab& S_fine = fine_level.get_new_data(State_Type);
#if CNS_USE_EB
    fine_level.flux_reg.Reflux(S_crse, *volfrac, S_fine, *fine_level.volfrac);
#else
    fine_level.flux_reg.Reflux(S_crse);
#endif

    if ((verbose != 0) && amrex::ParallelDescriptor::IOProcessor()) {
      amrex::Print() << " >> Reflux() at level " << level
                     << " : time = " << state[State_Type].curTime() << std::endl;
    }
  }

  if (level < parent->finestLevel()) {
    avgDown();
  }

  if (do_react && use_typical_vals_chem 
    && parent->levelSteps(0) % reset_typical_vals_int == 0) {
      set_typical_values_chem();
  }

  MultiFab& S = get_new_data(State_Type);
  MultiFab& IR = get_new_data(Reactions_Type);
  Real curtime = state[State_Type].curTime();
  Real dtlev = parent->dtLevel(level);
  Parm const* lparm = d_parm;
  ProbParm const* lprobparm = d_prob_parm;
  const auto geomdata = geom.data();
  auto const& sarrs = S.arrays();
  auto const& irarrs = IR.const_arrays();
  amrex::ParallelFor(S,
  [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept {
    // Do sth here, e.g. calculate time statistics
    prob_post_timestep(i, j, k, curtime, dtlev, sarrs[box_no], irarrs[box_no], 
                       geomdata, *lparm, *lprobparm);
  });
}

void
CNS::postCoarseTimeStep (Real time)
{ 
  // BL_PROFILE("CNS::postCoarseTimeStep()"); 
  // No need this because each function below has their own profile line

  AmrLevel::postCoarseTimeStep(time);

  printTotalandCheckNan(); //must do because this checks for nan as well
}

void
CNS::post_init (Real /*stop_time*/)
{
  // Initialise reaction
  if (do_react) { 
    if (use_typical_vals_chem) {
      set_typical_values_chem();
    }

    react_state(parent->cumTime(), parent->dtLevel(level), true);
  }
  
  if (level > 0) return;

  // Average data down from finer levels
  // so that conserved data is consistent between levels
  for (int k = parent->finestLevel()-1; k >= 0; --k) {
    getLevel(k).avgDown();
  }

  printTotalandCheckNan();
}

void
CNS::post_restart ()
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

  // Initialize reactor
  if (do_react) {
    if (chem_integrator == "ReactorNull") {
      Print() << "WARNING: turning on reactions while using ReactorNull. "
                 "Make sure this is intended." << std::endl;
    }

    reactor = pele::physics::reactions::ReactorBase::create(chem_integrator);
    reactor->init(1, 1);

    if (use_typical_vals_chem) {
      set_typical_values_chem();
    }
  }
  
#if (NUM_FIELD > 0)
  // Initialise random number
  Vector<IntVect> ref_ratio;
  for (int lv = 0; lv < parent->maxLevel(); ++lv) {
    ref_ratio.push_back(parent->refRatio(lv));
  }
  WienerProcess.init(AMREX_SPACEDIM, level, ref_ratio);
#endif

  MultiFab& S_new = get_new_data(State_Type);

  // Populate fields (when restarting from a different number of fields)
  if ((NUM_FIELD > 0) && do_restart_fields) {
    Print() << " >> Resetting stochastic fields state data ..." << std::endl;
    
    // Move aux variables
    MultiFab::Copy(S_new, S_new, NVAR, UFA, NUM_AUX, 0);

    // Copy mean to fields
    for (int nf = 1; nf < NUM_FIELD+1; ++nf) {
      MultiFab::Copy(S_new, S_new, 0, nf*NVAR, NVAR, 0);
    }

    if (do_react) {
      Print() << " >> Resetting stochastic fields reaction data ..." << std::endl;

      MultiFab& I_R = get_new_data(Reactions_Type);

      // Copy mean to fields
      for (int nf = 1; nf < NUM_FIELD+1; ++nf) {
        MultiFab::Copy(I_R, I_R, 0, nf*NREACT, NREACT, 0);
      }
    }
  }

  Parm const* lparm = d_parm;
  ProbParm const* lprobparm = d_prob_parm;
  const auto geomdata = geom.data();
  auto const& sarrs = S_new.arrays();
  amrex::ParallelFor(S_new,
  [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept {
    // Modify restarted data and/or add turbulence
    prob_post_restart(i, j, k, sarrs[box_no], geomdata, *lparm, *lprobparm);
  });
}

void
CNS::errorEst (TagBoxArray& tags, int /*clearval*/, int tagval, 
               Real time, int /*n_error_buf = 0*/, int /*ngrow*/)
{
  BL_PROFILE("CNS::errorEst()");

#if CNS_USE_EB
  if (refine_cutcells && (level < refine_cutcells_max_lev)) {
    const MultiFab& S_new = get_new_data(State_Type);
    TagCutCells(tags, S_new);
  }
#endif

  // for (int j = 0; j < errtagger.size(); ++j) {
  //     std::unique_ptr<MultiFab> mf;
  //     if (errtagger[j].Field() != std::string()) {
  //         mf = derive(errtagger[j].Field(), time, errtagger[j].NGrow());
  //     }
  //     errtagger[j](tags, mf.get(), clearval, tagval, time, level, geom);
  // }

  // Tag user specified box(es)
  if (!refine_boxes.empty()) {
    const int n_refine_boxes = refine_boxes.size();
    const auto problo = geom.ProbLoArray();
    const auto dx = geom.CellSizeArray();
    auto boxes = dp_refine_boxes;

    auto const& tagma = tags.arrays();
    ParallelFor(tags,
    [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept {
      RealVect pos {AMREX_D_DECL((i+0.5)*dx[0]+problo[0],
                                 (j+0.5)*dx[1]+problo[1],
                                 (k+0.5)*dx[2]+problo[2])};
      for (int irb = 0; irb < n_refine_boxes; ++irb) {
        if (boxes[irb].contains(pos) && (level < refine_boxes_max_lev[irb])) {
            tagma[box_no](i,j,k) = tagval;
        }
      }
    });
    // Gpu::streamSynchronize();
  }

  const MultiFab& S_new = get_new_data(State_Type);
  MultiFab Sborder(grids, dmap, NVAR, 1, MFInfo(), Factory());
  const Real cur_time = state[State_Type].curTime();
  FillPatch(*this, Sborder, Sborder.nGrow(), cur_time, State_Type, URHO, NVAR, 0);  
  
  auto const& tagma = tags.arrays();

  if (level < refine_dengrad_max_lev) {
    const Real threshold = refine_dengrad[level];

    MultiFab rho(Sborder, amrex::make_alias, URHO, 1);
    auto const& rhoma = rho.const_arrays();
    const auto dxinv = geom.InvCellSizeArray();
       
    ParallelFor(rho,
    [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept {
      cns_tag_grad(i, j, k, tagma[box_no], rhoma[box_no], dxinv, threshold, tagval);
    });
    // Gpu::streamSynchronize();
  }

  if (level < refine_velgrad_max_lev) {
    const Real threshold = refine_velgrad[level];

    for (MFIter mfi(tags, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
      const Box& bx = mfi.tilebox();
      
      FArrayBox velgrad(bx, 1);
      const int* bc = 0; //dummy
      cns_dervelgrad(bx, velgrad, 0, 1, Sborder[mfi], geom, cur_time, bc, level);
    
      auto tagarr = tags.array(mfi);
      auto mvarr = velgrad.const_array();
      ParallelFor(bx,
      [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        cns_tag_error(i, j, k, tagarr, mvarr, threshold, tagval);
      });
    }
    // Gpu::streamSynchronize();
  }

  if (level < refine_magvort_max_lev) {
    const Real threshold = refine_magvort[level];

    for (MFIter mfi(tags, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
      const Box& bx = mfi.tilebox();
      
      FArrayBox magvort(bx, 1);
      const int* bc = 0; //dummy
      cns_dermagvort(bx, magvort, 0, 1, Sborder[mfi], geom, cur_time, bc, level);
    
      auto tagarr = tags.array(mfi);
      auto mvarr = magvort.const_array();
      ParallelFor(bx,
      [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        cns_tag_error(i, j, k, tagarr, mvarr, threshold, tagval);
      });
    }
    // Gpu::streamSynchronize();
  }

#if (NUM_FIELD > 0) // cannot calculate tke if no SF
  if (level < refine_tke_max_lev) {
    const Real threshold = refine_tke[level];

    for (MFIter mfi(tags, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
      const Box& bx = mfi.tilebox();
      
      FArrayBox tke(bx, 1);
      const int* bc = 0; //dummy
      cns_dertke(bx, tke, 0, 1, S_new[mfi], geom, cur_time, bc, level);
    
      auto tagarr = tags.array(mfi);
      auto karr = tke.const_array();
      ParallelFor(bx,
      [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        cns_tag_error(i, j, k, tagarr, karr, threshold, tagval);
      });
    }
    // Gpu::streamSynchronize();
  }
#endif

  // // Tagging flame tracer
  // if (!flame_trac_name.empty()) {
  //     auto it = std::find(spec_names.begin(), spec_names.end(), flame_trac_name);

  //     if (it != spec_names.end()) {
  //         const auto idx = std::distance(spec_names.begin(), it);

  //         if (level < max_ftracer_lev) {
  //         const amrex::Real captured_ftracerr = tagging_parm->ftracerr;
  //         amrex::ParallelFor(
  //             tilebox, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
  //             tag_error(
  //                 i, j, k, tag_arr, S_derarr, captured_ftracerr, tagval);
  //             });
  //         }
  //         if (level < tagging_parm->max_ftracgrad_lev) {
  //         const amrex::Real captured_ftracgrad = tagging_parm->ftracgrad;
  //         amrex::ParallelFor(
  //             tilebox, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
  //             tag_graderror(
  //                 i, j, k, tag_arr, S_derarr, captured_ftracgrad, tagval);
  //             });
  //         }

  //     } else {
  //         amrex::Abort("Unknown species identified as flame_trac_name");
  //     }
  // }

  // Problem specific taggings
  Parm const* lparm = d_parm;
  ProbParm const* lprobparm = d_prob_parm;
  const auto geomdata = geom.data();
  auto tagarr = tags.arrays();
  auto const& sarrs = S_new.arrays();
  amrex::ParallelFor(tags,
  [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept {
    prob_tag_error(i, j, k, tagarr[box_no], sarrs[box_no], level, tagval, time, geomdata, *lparm, *lprobparm);
  });

  Gpu::streamSynchronize();
}

////////////////////////////////////////////////////////////////////////////
//    Helper functions that are in the CNS class but not in AMRLevel      //
////////////////////////////////////////////////////////////////////////////

/** 
 * \brief Print total and check for nan 
 */
void
CNS::printTotalandCheckNan () const
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
  Lazy::QueueReduction( [=] () mutable {
#endif
    ParallelDescriptor::ReduceRealSum(tot.data(), 5, ParallelDescriptor::IOProcessorNumber());
    if (verbose > 0) {
        amrex::Print().SetPrecision(17) << "\n[CNS] Total mass       = " << tot[0] << "\n"
                          AMREX_D_TERM( <<   "      Total x-momentum = " << tot[1] << "\n" ,
                                        <<   "      Total y-momentum = " << tot[2] << "\n" ,
                                        <<   "      Total z-momentum = " << tot[3] << "\n" )
                                        <<   "      Total energy     = " << tot[4] << "\n";
    }
#ifdef BL_LAZY
  });
#endif

  // Nan detector for soft exit
  if (isnan(tot[0]) || isnan(tot[4])) { signalStopJob = true; }
  amrex::ParallelDescriptor::ReduceBoolOr(signalStopJob);
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
int
CNS::okToContinue () 
{
  if (level > 0) {
    return 1;
  }

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

void
CNS::read_params ()
{
  ParmParse pp("cns");
  ParmParse pa("amr");

  pp.query("v", verbose);
  pp.query("cfl", cfl);
  pp.query("fixed_dt", fixed_dt);
  pp.query("dt_cutoff", dt_cutoff);

  pp.query("recon_char_var", recon_char_var);
  pp.query("char_sys", char_sys);
  if (char_sys != 0 && char_sys != 1) {
    amrex::Abort("CNS: char_sys must be 0 (speed of sound) or 1 (gamma)");
  }
  // recon_scheme
  //  1 -- piecewise constant
  //  2 -- piecewise linear
  //  3 -- WENO-Z 3rd order
  //  4 -- WENO-JS 5th order
  //  5 -- WENO-Z 5th order
  //  6 -- TENO 5
  pp.query("recon_scheme", recon_scheme); 
  // if (recon_scheme == 2) { // because we return to PLM near EBs, so give this option to users
  pp.query("limiter_theta", plm_theta); // MUSCL specific parameter
  // }

  Vector<int> lo_bc(AMREX_SPACEDIM), hi_bc(AMREX_SPACEDIM);
  pp.getarr("lo_bc", lo_bc, 0, AMREX_SPACEDIM);
  pp.getarr("hi_bc", hi_bc, 0, AMREX_SPACEDIM);
  for (int i = 0; i < AMREX_SPACEDIM; ++i) {
    phys_bc.setLo(i,lo_bc[i]);
    phys_bc.setHi(i,hi_bc[i]);
  }

  // pp.query("do_reflux", do_reflux); // How can you not do reflux?!

  pp.query("refine_cutcells", refine_cutcells);
  pp.query("refine_cutcells_max_lev", refine_cutcells_max_lev);

  pp.query("refine_dengrad_max_lev", refine_dengrad_max_lev);
  int numvals = pp.countval("refine_dengrad");
  int max_level; pa.query("max_level", max_level);
  if (numvals > 0) {
    if (numvals == max_level) {
      pp.queryarr("refine_dengrad", refine_dengrad);
    } else if (max_level > 0) {
      refine_dengrad.resize(max_level);
      Vector<Real> refine_tmp;
      pp.queryarr("refine_dengrad", refine_tmp);
      int lev = 0;
      for (; lev < std::min<int>(numvals, max_level); ++lev) {
        refine_dengrad[lev] = refine_tmp[lev];
      }
      for (; lev < max_level; ++lev) {
        refine_dengrad[lev] = refine_tmp[numvals-1];
      }
    }
  }

  pp.query("refine_velgrad_max_lev", refine_velgrad_max_lev);
  numvals = pp.countval("refine_velgrad");  
  if (numvals > 0) {
    if (numvals == max_level) {
      pp.queryarr("refine_velgrad", refine_velgrad);
    } else if (max_level > 0) {
      refine_velgrad.resize(max_level);
      Vector<Real> refine_tmp;
      pp.queryarr("refine_velgrad", refine_tmp);
      int lev = 0;
      for (; lev < std::min<int>(numvals, max_level); ++lev) {
        refine_velgrad[lev] = refine_tmp[lev];
      }
      for (; lev < max_level; ++lev) {
        refine_velgrad[lev] = refine_tmp[numvals-1];
      }
    }
  }

  pp.query("refine_magvort_max_lev", refine_magvort_max_lev);
  numvals = pp.countval("refine_magvort");
  if (numvals > 0) {
    if (numvals == max_level) {
      pp.queryarr("refine_magvort", refine_magvort);
    } else if (max_level > 0) {
      refine_magvort.resize(max_level);
      Vector<Real> refine_tmp;
      pp.queryarr("refine_magvort", refine_tmp);
      int lev = 0;
      for (; lev < std::min<int>(numvals, max_level); ++lev) {
        refine_magvort[lev] = refine_tmp[lev];
      }
      for (; lev < max_level; ++lev) {
        refine_magvort[lev] = refine_tmp[numvals-1];
      }
    } 
  }

#if (NUM_FIELD > 0) // cannot calculate tke if no SF
  pp.query("refine_tke_max_lev", refine_tke_max_lev);
  numvals = pp.countval("refine_tke");
  if (numvals > 0) {
    if (numvals == max_level) {
      pp.queryarr("refine_tke", refine_tke);
    } else if (max_level > 0) {
      refine_tke.resize(max_level);
      Vector<Real> refine_tmp;
      pp.queryarr("refine_tke", refine_tmp);
      int lev = 0;
      for (; lev < std::min<int>(numvals, max_level); ++lev) {
        refine_tke[lev] = refine_tmp[lev];
      }
      for (; lev < max_level; ++lev) {
        refine_tke[lev] = refine_tmp[numvals-1];
      }
    } 
  }
#endif

  int irefbox = 0;
  Vector<Real> refboxlo, refboxhi;
  int refbox_maxlev;
  while (pp.queryarr(("refine_box_lo_"+std::to_string(irefbox)).c_str(), refboxlo)) {
    pp.getarr(("refine_box_hi_"+std::to_string(irefbox)).c_str(), refboxhi);
    refine_boxes.emplace_back(refboxlo.data(), refboxhi.data());
    
    refbox_maxlev = 10;
    pp.query(("refine_box_max_level_"+std::to_string(irefbox)).c_str(), refbox_maxlev);
    refine_boxes_max_lev.emplace_back(refbox_maxlev);
    ++irefbox;
  }
  if (!refine_boxes.empty()) {
#ifdef AMREX_USE_GPU
    dp_refine_boxes = (RealBox*)The_Arena()->alloc(sizeof(RealBox)*refine_boxes.size());
    Gpu::htod_memcpy_async(dp_refine_boxes, refine_boxes.data(), sizeof(RealBox)*refine_boxes.size());
#else
    dp_refine_boxes = refine_boxes.data();
#endif
  }

  pp.query("do_visc", do_visc);
  pp.query("do_ext_src", do_ext_src);

  pp.query("do_react", do_react);
  pp.query("chem_integrator", chem_integrator);
  pp.query("rk_reaction_iter", rk_reaction_iter);
  pp.query("use_typical_vals_chem", use_typical_vals_chem);
  pp.query("reset_typical_vals_int", reset_typical_vals_int);
  pp.query("min_react_temp", min_react_temp);
  pp.query("clip_temp", clip_temp);
  pp.query("update_heat_release", update_heat_release);

  pp.query("do_restart_fields", do_restart_fields);
  pp.query("do_psgs", do_psgs);
  pp.query("do_vpdf", do_vpdf);
  pp.query("do_spdf", do_spdf);

  pp.query("do_les", do_les);
  pp.query("do_pasr", do_pasr);
  if (do_les) {
    if (AMREX_SPACEDIM == 1) amrex::Abort("CNS: LES not supporting 1D");
    pp.get("les_model", les_model_name);
    pp.query("Cs", Cs);
    pp.query("Pr_T", Pr_T);
    pp.query("Sc_T", Sc_T);
  }

#if CNS_USE_EB  
  // eb_weights_type:
  //   0 -- weights = 1             1 -- use_int_energy_as_eb_weights
  //   2 -- use_mass_as_eb_weights  3 -- use_volfrac_as_eb_weights
  pp.query("eb_weights_type", eb_weights_type);
  if (eb_weights_type < 0 || eb_weights_type > 3) {
    amrex::Abort("CNS: eb_weights_type must be 0, 1, 2 or 3");
  }

  pp.query("do_reredistribution", do_reredistribution);
  if (do_reredistribution != 0 && do_reredistribution != 1) {
    amrex::Abort("CNS: do_reredistibution must be 0 or 1");
  }

  pp.query("eb_no_slip", eb_no_slip);
  pp.query("eb_isothermal", eb_isothermal);
  if (eb_isothermal) {
    pp.get("eb_wall_temp", eb_wall_temp);
  }
#endif

  amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, h_parm, h_parm+1, d_parm);

  amrex::Gpu::streamSynchronize();
}

void
CNS::avgDown ()
{
  BL_PROFILE("CNS::avgDown()");

  if (level == parent->finestLevel()) return;

  auto& fine_lev = getLevel(level+1);
        MultiFab& S_crse =          get_new_data(State_Type);
  const MultiFab& S_fine = fine_lev.get_new_data(State_Type);

  MultiFab volume(S_fine.boxArray(), S_fine.DistributionMap(), 1, 0);
  geom.GetVolume(volume);

#if CNS_USE_EB
  amrex::EB_average_down(S_fine, S_crse, volume, fine_lev.volFrac(),
                         0, S_fine.nComp(), fine_ratio);
#else
  const amrex::Geometry& fgeom = getLevel(level + 1).geom;
  const amrex::Geometry& cgeom = geom;
  amrex::average_down(S_fine, S_crse, fgeom, cgeom, 
                      0, S_fine.nComp(), fine_ratio);
#endif

  if (do_react) {
          MultiFab& R_crse =          get_new_data(Reactions_Type);
    const MultiFab& R_fine = fine_lev.get_new_data(Reactions_Type);
#if CNS_USE_EB
    amrex::EB_average_down(R_fine, R_crse, volume, fine_lev.volFrac(),
                           0, R_fine.nComp(), fine_ratio);
#else
    amrex::average_down(R_fine, R_crse, fgeom, cgeom, 
                        0, R_fine.nComp(), fine_ratio);
#endif
  }
}

void
CNS::buildMetrics ()
{
  BL_PROFILE("CNS::buildMetrics()");
  
#if CNS_USE_EB    
  // make sure dx == dy == dz if use EB
  const Real* dx = geom.CellSize();

  if (AMREX_D_TERM(, std::abs(dx[0]-dx[1]) > 1.e-12*dx[0], 
                  || std::abs(dx[0]-dx[2]) > 1.e-12*dx[0])) {
    amrex::Abort("EB must have dx == dy == dz (for cut surface fluxes)\n");
  }

  if ((AMREX_SPACEDIM == 1) && (CNS_USE_EB)) {
    amrex::Abort("1D cannot do EB\n");
  }

  const auto& ebfactory = dynamic_cast<EBFArrayBoxFactory const&>(Factory());
  volfrac = &(ebfactory.getVolFrac());
  bndrycent = &(ebfactory.getBndryCent());
  areafrac = ebfactory.getAreaFrac();
  facecent = ebfactory.getFaceCent();

  Parm const* l_parm = d_parm;

  level_mask.clear();
  level_mask.define(grids,dmap,1,3);
  level_mask.BuildMask(geom.Domain(), geom.periodicity(),
                       l_parm->level_mask_covered,
                       l_parm->level_mask_notcovered,
                       l_parm->level_mask_physbnd,
                       l_parm->level_mask_interior);
#endif
}

Real
CNS::estTimeStep ()
{
  BL_PROFILE("CNS::estTimeStep()");

  if (fixed_dt > 0) { return fixed_dt; }

  const auto dx = geom.CellSizeArray();
  const MultiFab& S = get_new_data(State_Type);
  Parm const* lparm = d_parm;

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
  for (MFIter mfi(S,false); mfi.isValid(); ++mfi)
  {
    const Box& bx = mfi.tilebox();
    auto const& s_arr = S.array(mfi);
#if CNS_USE_EB
    const auto& flag = flags[mfi];
    if (flag.getType(bx) != FabType::covered)
#endif
    {
      reduce_op.eval(bx, reduce_data, [=]
      AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple {
        return cns_estdt_hydro(i, j, k, s_arr, dx, *lparm);
      });

      if (do_visc) {
        auto const* ltransparm = trans_parms.device_trans_parm();
        reduce_op.eval(bx, reduce_data, [=]
        AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple { 
          return cns_estdt_visc(i, j, k, s_arr, dx, *lparm, ltransparm);
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