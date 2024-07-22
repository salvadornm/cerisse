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
  store_in_checkpoint = do_react;
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
    cns_react_bcfill); // <--- need for all fields
  bndryfunc2.setRunOnGPU(true);

  desc_lst.setComponent(Reactions_Type, 0, react_name, react_bcs, bndryfunc2);

  // Setup Cost_Type
  desc_lst.addDescriptor(Cost_Type, IndexType::TheCellType(), StateDescriptor::Point,
                         0, 1, &pc_interp);
  desc_lst.setComponent(Cost_Type, 0, "Cost", bc, bndryfunc);

  // SNM: commeted
  // assert(num_state_data_types == desc_lst.size());

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

  //
  // Dynamically generated error tagging functions
  //
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
      const std::string field = "mag_vort";
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
