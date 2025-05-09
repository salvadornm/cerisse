#include <AMReX_Derive.H>

#include "CNS.H"
#include "derive.H"
#include "prob.H"

using namespace amrex;

ProbParm* CNS::h_prob_parm = nullptr;
ProbParm* CNS::d_prob_parm = nullptr;
#ifdef USE_PROB_PARM_HOST
ProbParmHost* CNS::prob_parm_host = nullptr;
#endif

using BndryFunc = StateDescriptor::BndryFunc;

// Components are:
// Interior, Inflow = UserBC, Outflow, Symmetry = SlipWall, SlipWall, NoSlipWall,
// UserBC
static int scalar_bc[] = {
  BCType::int_dir,      BCType::ext_dir,      BCType::foextrap, BCType::reflect_even,
  BCType::reflect_even, BCType::reflect_even, BCType::ext_dir};

static int norm_vel_bc[] = {
  BCType::int_dir,     BCType::ext_dir,     BCType::foextrap, BCType::reflect_odd,
  BCType::reflect_odd, BCType::reflect_odd, BCType::ext_dir};

static int tang_vel_bc[] = {
  BCType::int_dir,      BCType::ext_dir,     BCType::foextrap, BCType::reflect_even,
  BCType::reflect_even, BCType::reflect_odd, BCType::ext_dir};

static int react_src_bc[] = {BCType::int_dir,      BCType::reflect_even,
                             BCType::reflect_even, BCType::reflect_even,
                             BCType::reflect_even, BCType::reflect_even,
                             BCType::reflect_even};

static void set_scalar_bc(BCRec& bc, const BCRec& phys_bc)
{
  const int* lo_bc = phys_bc.lo();
  const int* hi_bc = phys_bc.hi();
  for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
    bc.setLo(dir, scalar_bc[lo_bc[dir]]);
    bc.setHi(dir, scalar_bc[hi_bc[dir]]);
  }
}

static void set_x_vel_bc(BCRec& bc, const BCRec& phys_bc)
{
  const int* lo_bc = phys_bc.lo();
  const int* hi_bc = phys_bc.hi();
  AMREX_D_TERM(
    bc.setLo(0, norm_vel_bc[lo_bc[0]]); bc.setHi(0, norm_vel_bc[hi_bc[0]]);
    , bc.setLo(1, tang_vel_bc[lo_bc[1]]); bc.setHi(1, tang_vel_bc[hi_bc[1]]);
    , bc.setLo(2, tang_vel_bc[lo_bc[2]]); bc.setHi(2, tang_vel_bc[hi_bc[2]]);)
}

static void set_y_vel_bc(BCRec& bc, const BCRec& phys_bc)
{
  const int* lo_bc = phys_bc.lo();
  const int* hi_bc = phys_bc.hi();
  AMREX_D_TERM(
    bc.setLo(0, tang_vel_bc[lo_bc[0]]); bc.setHi(0, tang_vel_bc[hi_bc[0]]);
    , bc.setLo(1, norm_vel_bc[lo_bc[1]]); bc.setHi(1, norm_vel_bc[hi_bc[1]]);
    , bc.setLo(2, tang_vel_bc[lo_bc[2]]); bc.setHi(2, tang_vel_bc[hi_bc[2]]);)
}

static void set_z_vel_bc(BCRec& bc, const BCRec& phys_bc)
{
  const int* lo_bc = phys_bc.lo();
  const int* hi_bc = phys_bc.hi();
  AMREX_D_TERM(
    bc.setLo(0, tang_vel_bc[lo_bc[0]]); bc.setHi(0, tang_vel_bc[hi_bc[0]]);
    , bc.setLo(1, tang_vel_bc[lo_bc[1]]); bc.setHi(1, tang_vel_bc[hi_bc[1]]);
    , bc.setLo(2, norm_vel_bc[lo_bc[2]]); bc.setHi(2, norm_vel_bc[hi_bc[2]]);)
}

static void set_react_src_bc(BCRec& bc, const BCRec& phys_bc)
{
  const int* lo_bc = phys_bc.lo();
  const int* hi_bc = phys_bc.hi();
  for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
    bc.setLo(dir, react_src_bc[lo_bc[dir]]);
    bc.setHi(dir, react_src_bc[hi_bc[dir]]);
  }
}


void CNS::read_params()
{
  ParmParse pp("cns");
  ParmParse pa("amr");

  pp.query("v", verbose);
  pp.query("cfl", cfl);
  pp.query("fixed_dt", fixed_dt);
  pp.query("dt_cutoff", dt_cutoff);
  pp.query("dt_max_change", dt_max_change);
  pp.query("check_message_int", check_message_int);

#if CNS_USE_EB
  // eb_weight in redist
  pp.query("eb_weight", eb_weight);

  pp.query("redistribution_type", redistribution_type);
  if (redistribution_type != "StateRedist" && redistribution_type != "FluxRedist" &&
      redistribution_type != "NoRedist") {
    amrex::Abort(
      "CNS: redistribution_type must be StateRedist or FluxRedist or NoRedist");
  }

  pp.query("eb_no_slip", eb_no_slip);
  pp.query("eb_isothermal", eb_isothermal);
  if (eb_isothermal) { pp.get("eb_wall_temp", eb_wall_temp); }
  pp.query("eb_recon_mode", eb_recon_mode); // 0: fill stencil, 1: return to PLM

  pp.query("eb_wall_model", eb_wall_model);
#endif

  pp.query("recon_char_var", recon_char_var);
  pp.query("char_sys", char_sys);
  if (char_sys != 0 && char_sys != 1) {
    amrex::Abort("CNS: char_sys must be 0 (speed of sound) or 1 (gamma)");
  }
  // recon_scheme
  //  1 -- piecewise constant; 2 -- piecewise linear; 3 -- WENO-Z 3rd order
  //  4 -- WENO-JS 5th order;  5 -- WENO-Z 5th order; 6 -- TENO 5
  pp.query("recon_scheme", recon_scheme);
  if (recon_scheme == 2
#if CNS_USE_EB
      // because return to PLM near EBs, so give this option to users
      || eb_recon_mode == 1
#endif
  ) {
    pp.query("limiter_theta", plm_theta); // MUSCL limiter parameter
  }
  pp.query("use_hybrid_scheme", use_hybrid_scheme);
  pp.query("teno_cutoff", teno_cutoff);

  Vector<int> lo_bc(AMREX_SPACEDIM), hi_bc(AMREX_SPACEDIM);
  pp.getarr("lo_bc", lo_bc, 0, AMREX_SPACEDIM);
  pp.getarr("hi_bc", hi_bc, 0, AMREX_SPACEDIM);
  for (int i = 0; i < AMREX_SPACEDIM; ++i) {
    phys_bc.setLo(i, lo_bc[i]);
    phys_bc.setHi(i, hi_bc[i]);
  }

  pp.query("do_reflux", do_reflux); // How can you not do reflux?!

  pp.query("refine_dengrad_max_lev", refine_dengrad_max_lev);
  int numvals = pp.countval("refine_dengrad");
  int max_level;
  pa.query("max_level", max_level);
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
        refine_dengrad[lev] = refine_tmp[numvals - 1];
      }
    }
  }

#if CNS_USE_EB
  pp.query("refine_cutcells", refine_cutcells);
  pp.query("refine_cutcells_max_lev", refine_cutcells_max_lev);
#endif

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
        refine_velgrad[lev] = refine_tmp[numvals - 1];
      }
    }
  }

  pp.query("refine_presgrad_max_lev", refine_presgrad_max_lev);
  numvals = pp.countval("refine_presgrad");
  if (numvals > 0) {
    if (numvals == max_level) {
      pp.queryarr("refine_presgrad", refine_presgrad);
    } else if (max_level > 0) {
      refine_presgrad.resize(max_level);
      Vector<Real> refine_tmp;
      pp.queryarr("refine_presgrad", refine_tmp);
      int lev = 0;
      for (; lev < std::min<int>(numvals, max_level); ++lev) {
        refine_presgrad[lev] = refine_tmp[lev];
      }
      for (; lev < max_level; ++lev) {
        refine_presgrad[lev] = refine_tmp[numvals - 1];
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
        refine_magvort[lev] = refine_tmp[numvals - 1];
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
      for (; lev < max_level; ++lev) { refine_tke[lev] = refine_tmp[numvals - 1]; }
    }
  }
#endif

  int irefbox = 0;
  Vector<Real> refboxlo, refboxhi;
  int refbox_maxlev;
  while (
    pp.queryarr(("refine_box_lo_" + std::to_string(irefbox)).c_str(), refboxlo)) {
    pp.getarr(("refine_box_hi_" + std::to_string(irefbox)).c_str(), refboxhi);
    refine_boxes.emplace_back(refboxlo.data(), refboxhi.data());

    refbox_maxlev = 10;
    pp.query(("refine_box_max_level_" + std::to_string(irefbox)).c_str(),
             refbox_maxlev);
    refine_boxes_max_lev.emplace_back(refbox_maxlev);
    ++irefbox;
  }
  if (!refine_boxes.empty()) {
#ifdef AMREX_USE_GPU
    dp_refine_boxes =
      (RealBox*)The_Arena()->alloc(sizeof(RealBox) * refine_boxes.size());
    Gpu::htod_memcpy_async(dp_refine_boxes, refine_boxes.data(),
                           sizeof(RealBox) * refine_boxes.size());
#else
    dp_refine_boxes = refine_boxes.data();
#endif
  }

  if (pp.queryarr("buffer_box_lo", refboxlo)) {
    pp.getarr("buffer_box_hi", refboxhi);
    buffer_box.setLo(refboxlo.data());
    buffer_box.setHi(refboxhi.data());
    if (!buffer_box.ok()) { amrex::Abort("CNS: buffer_box has negative volume"); }
  }

  pp.query("do_hydro", do_hydro);
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
  pp.query("do_pd_model", do_pd_model);
  pp.query("do_vpdf", do_vpdf);
  pp.query("do_spdf", do_spdf);
  pp.query("do_species_langevin", do_species_langevin);

  pp.query("amr_interp_order", amr_interp_order);
  if (amr_interp_order == 3 && amrex::SpaceDim == 1) {
    amrex::Abort("CNS: quadratic_interp does not support 1D.");
  }
  pp.query("rk_order", rk_order); // 1-4
  if (rk_order != 1 && rk_order != 2) {
    if (do_react || do_psgs || do_pd_model || do_vpdf || do_spdf) {
      amrex::Abort("CNS: rk_order must be 1 or 2 for reaction and/or SF");
    }
  }

  pp.query("do_les", do_les);
  pp.query("do_pasr", do_pasr);
  if (do_les || do_pasr) {
    if (AMREX_SPACEDIM == 1) amrex::Abort("CNS: LES not supporting 1D");
    pp.query("C_s", Cs);
    pp.query("C_I", C_I);
    pp.query("Pr_T", Pr_T);
    pp.query("Sc_T", Sc_T);
    pp.query("Cm", Cm);
  }

  pp.query("do_nscbc", do_nscbc);
  pp.query("nscbc_relax_p", nscbc_relax_p);
  pp.query("nscbc_relax_u", nscbc_relax_u);
  pp.query("nscbc_relax_T", nscbc_relax_T);
  pp.query("ambient_p", ambient_p);
}

void CNS::read_errtags() {
  ParmParse pp("cns");
  Vector<std::string> refinement_indicators;
  pp.queryarr("refinement_indicators", refinement_indicators);

  for (int i = 0; i < refinement_indicators.size(); ++i) {
    std::string ref_prefix = "cns." + refinement_indicators[i];
    ParmParse ppr(ref_prefix);

    RealBox realbox;
    if (ppr.countval("in_box_lo")) {
      std::vector<Real> box_lo(amrex::SpaceDim), box_hi(amrex::SpaceDim);
      ppr.getarr("in_box_lo", box_lo, 0, box_lo.size());
      ppr.getarr("in_box_hi", box_hi, 0, box_hi.size());
      realbox = RealBox(&(box_lo[0]), &(box_hi[0]));
    }

    AMRErrorTagInfo info;
    if (realbox.ok()) { info.SetRealBox(realbox); }
    if (ppr.countval("start_time") > 0) {
      Real min_time;
      ppr.get("start_time", min_time);
      info.SetMinTime(min_time);
    }
    if (ppr.countval("end_time") > 0) {
      Real max_time;
      ppr.get("end_time", max_time);
      info.SetMaxTime(max_time);
    }
    if (ppr.countval("max_level") > 0) {
      int max_level;
      ppr.get("max_level", max_level);
      info.SetMaxLevel(max_level);
    }

    if (ppr.countval("value_greater")) {
      Real value;
      ppr.get("value_greater", value);
      std::string field;
      ppr.get("field_name", field);
      errtags.push_back(AMRErrorTag(value, AMRErrorTag::GREATER, field, info));
    } else if (ppr.countval("value_less")) {
      Real value;
      ppr.get("value_less", value);
      std::string field;
      ppr.get("field_name", field);
      errtags.push_back(AMRErrorTag(value, AMRErrorTag::LESS, field, info));
    } else if (ppr.countval("vorticity_greater")) {
      Real value;
      ppr.get("vorticity_greater", value);
      const std::string field = "magvort";
      errtags.push_back(AMRErrorTag(value, AMRErrorTag::VORT, field, info));
    } else if (ppr.countval("adjacent_difference_greater")) {
      Real value;
      ppr.get("adjacent_difference_greater", value);
      std::string field;
      ppr.get("field_name", field);
      errtags.push_back(AMRErrorTag(value, AMRErrorTag::GRAD, field, info));
    } else if (realbox.ok()) {
      errtags.push_back(AMRErrorTag(info));
    } else {
      Abort(std::string("Unrecognised refinement indicator for " +
                        refinement_indicators[i])
              .c_str());
    }
  }
}

void CNS::variableSetUp()
{
  h_prob_parm = new ProbParm{}; // This lives on host
#ifdef USE_PROB_PARM_HOST
  prob_parm_host = new ProbParmHost{}; // This lives on host
#endif
  d_prob_parm =
    (ProbParm*)The_Arena()->alloc(sizeof(ProbParm)); // This lives on device

  trans_parms.allocate();              // PelePhysics trans_parms
  turb_inflow.init(DefaultGeometry()); // PelePhysics turb_inflow

  read_params();

  // - pc_interp       : Piecewise constant interpolation.
  // - cell_bilinear_interp: Linear interpolation for cell-centered data. No 3D.
  // - quadratic_interp: Quadratic polynomial interpolation. No 1D.
  // - lincc_interp    : Dimension-by-dimension linear interpolation with MC limiter.
  //   (PeleC default)   For multi-component data, the strictest limiter is used for
  //                     all components. Conservative for both Cartesian and
  //                     curvilinear coordinates.
  // - cell_cons_interp: Similar to lincc_interp, but 1) the interpolations for each
  //   (CNS default)     component are independent of each other; 2) after the
  //                     dimension-by-dimension linear interpolation with limiting,
  //                     there is a further limiting to ensure no new min or max are
  //                     created in fine cells.
  // - protected_interp: Similar to lincc_interp. Additionally, it has a  function
  //                     one can call to ensure no values are negative.
  // - quartic_interp  : Quartic polynomial conservative interpolation.
  Interpolater* interp;
#if CNS_USE_EB
    // EB only supports lincc or cell_cons
    interp = &eb_lincc_interp;
    // interp = &eb_cell_cons_interp;
#else
  if (amr_interp_order == 1) {
    interp = &pc_interp;
  } else if (amr_interp_order == 3) {
    interp = &quadratic_interp;
  } else if (amr_interp_order == 4) {
    interp = &quartic_interp;
  } else {
    // default: second-order limited interp
    interp = &lincc_interp;
    // interp = &cell_cons_interp;
  }
#endif

  // Setup State_Type
  bool state_data_extrap = false;
  bool store_in_checkpoint = true;
  desc_lst.addDescriptor(State_Type, IndexType::TheCellType(),
                         StateDescriptor::Point, NUM_GROW, LEN_STATE, interp,
                         state_data_extrap, store_in_checkpoint);

  Vector<amrex::BCRec> bcs(LEN_STATE);
  Vector<std::string> name(LEN_STATE);
  BCRec bc;

  // Get the species names from the network model
  pele::physics::eos::speciesNames<pele::physics::eos::Fuego>(spec_names);

  if (ParallelDescriptor::IOProcessor()) {
    amrex::Print() << NUM_REACTIONS << " Reactions in mechanism \n";
    // Print species names
    amrex::Print() << NUM_SPECIES << " Species: " << std::endl;
    for (int i = 0; i < NUM_SPECIES; i++) {
      amrex::Print() << spec_names[i] << ' ' << ' ';
    }
    amrex::Print() << std::endl;

    // Print number of fields and aux
    amrex::Print() << NUM_FIELD << " Fields, " << NUM_AUX << " Auxiliary Variables, "
                   << LEN_STATE << LEN_REACT << LEN_PRIM << LEN_COEF << '\n';
  }

  int cnt = 0; // variable counter

  // Mean field
  set_scalar_bc(bc, phys_bc);
  bcs[cnt] = bc;
  name[cnt] = "density";
  cnt++;
  set_x_vel_bc(bc, phys_bc);
  bcs[cnt] = bc;
  name[cnt] = "xmom";
  cnt++;
  set_y_vel_bc(bc, phys_bc);
  bcs[cnt] = bc;
  name[cnt] = "ymom";
  cnt++;
  set_z_vel_bc(bc, phys_bc);
  bcs[cnt] = bc;
  name[cnt] = "zmom";
  cnt++;
  set_scalar_bc(bc, phys_bc);
  bcs[cnt] = bc;
  name[cnt] = "rho_E";
  cnt++;
  set_scalar_bc(bc, phys_bc);
  bcs[cnt] = bc;
  name[cnt] = "Temp";
  cnt++;
  for (int i = 0; i < NUM_SPECIES; ++i) {
    set_scalar_bc(bc, phys_bc);
    bcs[cnt] = bc;
    name[cnt] = "rho_" + spec_names[i];
    cnt++;
  }

  // For each field
  for (int nf = 0; nf < NUM_FIELD; ++nf) {
    set_scalar_bc(bc, phys_bc);
    bcs[cnt] = bc;
    name[cnt] = "density_Field" + std::to_string(nf);
    cnt++;
    set_x_vel_bc(bc, phys_bc);
    bcs[cnt] = bc;
    name[cnt] = "xmom_Field" + std::to_string(nf);
    cnt++;
    set_y_vel_bc(bc, phys_bc);
    bcs[cnt] = bc;
    name[cnt] = "ymom_Field" + std::to_string(nf);
    cnt++;
    set_z_vel_bc(bc, phys_bc);
    bcs[cnt] = bc;
    name[cnt] = "zmom_Field" + std::to_string(nf);
    cnt++;
    set_scalar_bc(bc, phys_bc);
    bcs[cnt] = bc;
    name[cnt] = "rho_E_Field" + std::to_string(nf);
    cnt++;
    set_scalar_bc(bc, phys_bc);
    bcs[cnt] = bc;
    name[cnt] = "Temp_Field" + std::to_string(nf);
    cnt++;
    for (int i = 0; i < NUM_SPECIES; ++i) {
      set_scalar_bc(bc, phys_bc);
      bcs[cnt] = bc;
      name[cnt] = "rho_" + spec_names[i] + "_Field" + std::to_string(nf);
      cnt++;
    }
  }

#if NUM_AUX > 0
  // Get AUX names
  amrex::Vector<std::string> aux_name;
  prob_get_aux_name(aux_name);
  for (int i = 0; i < NUM_AUX; ++i) {
    set_scalar_bc(bc, phys_bc);
    bcs[cnt] = bc;
    if (aux_name.size() == NUM_AUX) {
      name[cnt] = aux_name[i];
    } else {
      name[cnt] = "aux_" + std::to_string(i);
    }
    cnt++;
  }
#endif

  StateDescriptor::BndryFunc bndryfunc(cns_bcfill);
  bndryfunc.setRunOnGPU(true);

  desc_lst.setComponent(State_Type, 0, name, bcs, bndryfunc);

  // Setup React_Type
  store_in_checkpoint = false;
  desc_lst.addDescriptor(Reactions_Type, IndexType::TheCellType(),
                         StateDescriptor::Point, 0, LEN_REACT, interp,
                         state_data_extrap, store_in_checkpoint);

  Vector<amrex::BCRec> react_bcs(LEN_REACT);
  Vector<std::string> react_name(LEN_REACT);

  cnt = 0;

  // Mean field
  for (int i = 0; i < NUM_SPECIES; ++i) {
    set_react_src_bc(bc, phys_bc);
    react_bcs[cnt] = bc;
    react_name[cnt] = "rho_omega_" + spec_names[i];
    cnt++;
  }
  set_react_src_bc(bc, phys_bc);
  react_bcs[cnt] = bc;
  react_name[cnt] = "heatRelease";
  cnt++;

  // For each field
  for (int nf = 0; nf < NUM_FIELD; ++nf) {
    for (int i = 0; i < NUM_SPECIES; ++i) {
      set_react_src_bc(bc, phys_bc);
      react_bcs[cnt] = bc;
      react_name[cnt] = "rho_omega_" + spec_names[i] + "_Field" + std::to_string(nf);
      cnt++;
    }
    set_react_src_bc(bc, phys_bc);
    react_bcs[cnt] = bc;
    react_name[cnt] = "heatRelease_Field" + std::to_string(nf);
    cnt++;
  }

  StateDescriptor::BndryFunc bndryfunc2(
    cns_null_bcfill); // <--- need for all fields
  bndryfunc2.setRunOnGPU(true);

  desc_lst.setComponent(Reactions_Type, 0, react_name, react_bcs, bndryfunc2);

  // Setup Cost_Type
  desc_lst.addDescriptor(Cost_Type, IndexType::TheCellType(), StateDescriptor::Point,
                         0, 1, &pc_interp, state_data_extrap, store_in_checkpoint);
  desc_lst.setComponent(Cost_Type, 0, "Cost", bc, bndryfunc2);

  Print() << desc_lst.size() << " Data Types:" << std::endl;
  for (int typ = 0; typ < desc_lst.size(); typ++) {
    Print() << typ << " - ";
    const StateDescriptor& desc = desc_lst[typ];
    for (int n = 0; n < desc.nComp(); n++) { Print() << desc.name(n) << " "; }
    Print() << "(" << desc.nComp() << ")" << std::endl;
  }

  StateDescriptor::setBndryFuncThreadSafety(true);

  // DEFINE DERIVED QUANTITIES (Derive from MEAN field)
  // Pressure
  derive_lst.add("pressure", IndexType::TheCellType(), 1, cns_derpres,
                 DeriveRec::TheSameBox);
  derive_lst.addComponent("pressure", desc_lst, State_Type, 0, NVAR);

  // Internal energy
  derive_lst.add("eint", IndexType::TheCellType(), 1, cns_dereint,
                 DeriveRec::TheSameBox);
  derive_lst.addComponent("eint", desc_lst, State_Type, 0, NVAR);

  // Velocities
  derive_lst.add("velocity", IndexType::TheCellType(), amrex::SpaceDim,
                 {AMREX_D_DECL("x_velocity", "y_velocity", "z_velocity")},
                 cns_dervel, DeriveRec::TheSameBox);
  derive_lst.addComponent("velocity", desc_lst, State_Type, 0, NVAR);

  // Mach number
  derive_lst.add("MachNumber", IndexType::TheCellType(), 1, cns_dermachnumber,
                 DeriveRec::TheSameBox);
  derive_lst.addComponent("MachNumber", desc_lst, State_Type, 0, NVAR);

  // Note: All GrowBoxByOne derived quantities need all NVAR states because bcnormal
  // may need them (e.g. isothermal wall)

  // Vorticity
  derive_lst.add("magvort", IndexType::TheCellType(), 1, cns_dermagvort,
                 DeriveRec::GrowBoxByOne);
  derive_lst.addComponent("magvort", desc_lst, State_Type, 0, NVAR);
  // Numerical schlieren
  derive_lst.add("divu", IndexType::TheCellType(), 1, cns_derdivu,
                 DeriveRec::GrowBoxByOne);
  derive_lst.addComponent("divu", desc_lst, State_Type, 0, NVAR);

  derive_lst.add("divrho", IndexType::TheCellType(), 1, cns_derdivrho,
                 DeriveRec::GrowBoxByOne);
  derive_lst.addComponent("divrho", desc_lst, State_Type, 0, NVAR);

  derive_lst.add("shock_sensor", IndexType::TheCellType(), 1, cns_dershocksensor,
                 DeriveRec::GrowBoxByOne);
  derive_lst.addComponent("shock_sensor", desc_lst, State_Type, 0, NVAR);

  // External source term
  // derive_lst.add("ext_src", IndexType::TheCellType(), 1,
  //                CNS::cns_derextsrc, DeriveRec::TheSameBox);
  // derive_lst.addComponent("ext_src", desc_lst, State_Type, URHO, NVAR);

  // Cp and Cv
  derive_lst.add("cp", IndexType::TheCellType(), 1, cns_dercp,
                 DeriveRec::TheSameBox);
  derive_lst.addComponent("cp", desc_lst, State_Type, 0, NVAR);

  derive_lst.add("cv", IndexType::TheCellType(), 1, cns_dercv,
                 DeriveRec::TheSameBox);
  derive_lst.addComponent("cv", desc_lst, State_Type, 0, NVAR);

  // Transport coefficients
  Vector<std::string> trans_coef_names(NUM_SPECIES + 3);
  for (int i = 0; i < NUM_SPECIES; i++) {
    trans_coef_names[i] = "D_" + spec_names[i];
  }
  trans_coef_names[NUM_SPECIES] = "viscosity";
  trans_coef_names[NUM_SPECIES + 1] = "bulk_viscosity";
  trans_coef_names[NUM_SPECIES + 2] = "conductivity";

  derive_lst.add("transport_coef", IndexType::TheCellType(), NUM_SPECIES + 3,
                 trans_coef_names, CNS::cns_dertranscoef, DeriveRec::TheSameBox);
  derive_lst.addComponent("transport_coef", desc_lst, State_Type, 0, NVAR);

  // Derived mass and mole fractions (for all fields)
  amrex::Vector<std::string> massfrac_names((NUM_FIELD + 1) * NUM_SPECIES);
  amrex::Vector<std::string> molefrac_names((NUM_FIELD + 1) * NUM_SPECIES);
  for (int i = 0; i < NUM_SPECIES; i++) {
    massfrac_names[i] = "Y_" + spec_names[i];
    molefrac_names[i] = "X_" + spec_names[i];
  }
  derive_lst.add("massfrac", IndexType::TheCellType(), (NUM_FIELD + 1) * NUM_SPECIES,
                 massfrac_names, cns_dermassfrac, DeriveRec::TheSameBox);
  derive_lst.addComponent("massfrac", desc_lst, State_Type, 0, LEN_STATE);
  derive_lst.add("molefrac", IndexType::TheCellType(), (NUM_FIELD + 1) * NUM_SPECIES,
                 molefrac_names, cns_dermolefrac, DeriveRec::TheSameBox);
  derive_lst.addComponent("molefrac", desc_lst, State_Type, 0, LEN_STATE);

  // LES model
  if (do_les) {
    Vector<std::string> turb_visc_names = {"mu_t", "xi_t"};
    derive_lst.add("turb_viscosity", IndexType::TheCellType(), 2, turb_visc_names,
                   cns_derturbvisc, DeriveRec::GrowBoxByOne);
    derive_lst.addComponent("turb_viscosity", desc_lst, State_Type, URHO, NVAR);
  }

#if NUM_FIELD > 0
  // Derived variance of velocity, species, pressure
  amrex::Vector<std::string> rs_names = {"R11", "R12", "R22", "R13", "R23", "R33"};
  derive_lst.add("reynolds_stress", IndexType::TheCellType(), AMREX_D_PICK(1, 3, 6),
                 rs_names, cns_dervaru, DeriveRec::TheSameBox);
  derive_lst.addComponent("reynolds_stress", desc_lst, State_Type, 0, LEN_STATE);

  amrex::Vector<std::string> vary_names(NUM_SPECIES);
  for (int i = 0; i < NUM_SPECIES; i++) { vary_names[i] = "var_Y_" + spec_names[i]; }
  derive_lst.add("var_Y", IndexType::TheCellType(), NUM_SPECIES, vary_names,
                 cns_dervary, DeriveRec::TheSameBox);
  derive_lst.addComponent("var_Y", desc_lst, State_Type, 0, LEN_STATE);

  derive_lst.add("var_p", IndexType::TheCellType(), 1, cns_dervarp,
                 DeriveRec::TheSameBox);
  derive_lst.addComponent("var_p", desc_lst, State_Type, 0, LEN_STATE);

  derive_lst.add("var_T", IndexType::TheCellType(), 1, cns_dervarT,
                 DeriveRec::TheSameBox);
  derive_lst.addComponent("var_T", desc_lst, State_Type, 0, LEN_STATE);
#endif

  // Dynamically generated error tagging functions
  CNS::read_errtags();
}

void CNS::variableCleanUp()
{
  delete h_prob_parm;
  The_Arena()->free(d_prob_parm);
#ifdef USE_PROB_PARM_HOST
  delete prob_parm_host;
#endif
  desc_lst.clear();
  derive_lst.clear();

  trans_parms.deallocate();
  pmf_data.deallocate();

#ifdef AMREX_USE_GPU
  The_Arena()->free(dp_refine_boxes);
#endif
}
