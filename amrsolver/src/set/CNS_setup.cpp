#include <CNS.h>
#include <CNS_K.h>
#include <prob.h>

using namespace amrex;

int CNS::num_state_data_types = 0;

static Box the_same_box(const Box& b) { return b; }
// static Box grow_box_by_one (const Box& b) { return amrex::grow(b,1); }

using BndryFunc = StateDescriptor::BndryFunc;

//
// Components are:
//  Interior, Inflow, Outflow, Symmetry, SlipWall, NoSlipWall, User defined
//
static int scalar_bc[] = {BCType::int_dir,      BCType::ext_dir,
                          BCType::foextrap,     BCType::reflect_even,
                          BCType::reflect_even, BCType::reflect_even,
                          BCType::ext_dir};

static int norm_vel_bc[] = {BCType::int_dir,     BCType::ext_dir,
                            BCType::foextrap,    BCType::reflect_odd,
                            BCType::reflect_odd, BCType::reflect_odd,
                            BCType::ext_dir};

static int tang_vel_bc[] = {BCType::int_dir,      BCType::ext_dir,
                            BCType::foextrap,     BCType::reflect_even,
                            BCType::reflect_even, BCType::reflect_odd,
                            BCType::ext_dir};

static void set_scalar_bc(BCRec& bc, const BCRec* phys_bc) {
  const int* lo_bc = phys_bc->lo();
  const int* hi_bc = phys_bc->hi();
  for (int i = 0; i < AMREX_SPACEDIM; i++) {
    bc.setLo(i, scalar_bc[lo_bc[i]]);
    bc.setHi(i, scalar_bc[hi_bc[i]]);
  }
}

static void set_x_vel_bc(BCRec& bc, const BCRec* phys_bc) {
  const int* lo_bc = phys_bc->lo();
  const int* hi_bc = phys_bc->hi();
  bc.setLo(0, norm_vel_bc[lo_bc[0]]);
  bc.setHi(0, norm_vel_bc[hi_bc[0]]);
#if (AMREX_SPACEDIM >= 2)
  bc.setLo(1, tang_vel_bc[lo_bc[1]]);
  bc.setHi(1, tang_vel_bc[hi_bc[1]]);
#endif
#if (AMREX_SPACEDIM == 3)
  bc.setLo(2, tang_vel_bc[lo_bc[2]]);
  bc.setHi(2, tang_vel_bc[hi_bc[2]]);
#endif
}

static void set_y_vel_bc(BCRec& bc, const BCRec* phys_bc) {
  const int* lo_bc = phys_bc->lo();
  const int* hi_bc = phys_bc->hi();
  bc.setLo(0, tang_vel_bc[lo_bc[0]]);
  bc.setHi(0, tang_vel_bc[hi_bc[0]]);
#if (AMREX_SPACEDIM >= 2)
  bc.setLo(1, norm_vel_bc[lo_bc[1]]);
  bc.setHi(1, norm_vel_bc[hi_bc[1]]);
#endif
#if (AMREX_SPACEDIM == 3)
  bc.setLo(2, tang_vel_bc[lo_bc[2]]);
  bc.setHi(2, tang_vel_bc[hi_bc[2]]);
#endif
}

#if (AMREX_SPACEDIM == 3)
static void set_z_vel_bc(BCRec& bc, const BCRec* phys_bc) {
  const int* lo_bc = phys_bc->lo();
  const int* hi_bc = phys_bc->hi();
  bc.setLo(0, tang_vel_bc[lo_bc[0]]);
  bc.setHi(0, tang_vel_bc[hi_bc[0]]);
  bc.setLo(1, tang_vel_bc[lo_bc[1]]);
  bc.setHi(1, tang_vel_bc[hi_bc[1]]);
  bc.setLo(2, norm_vel_bc[lo_bc[2]]);
  bc.setHi(2, norm_vel_bc[hi_bc[2]]);
}
#endif

void CNS::variableSetUp() {
  // Closures and Problem structures (available on both CPU and GPU)
  CNS::h_prob_closures = new PROB::ProbClosures{};
  CNS::h_prob_parm = new PROB::ProbParm{};
  CNS::h_phys_bc = new BCRec{};
#ifdef AMREX_USE_GPU
  CNS::d_prob_closures =
      (PROB::ProbClosures*)The_Arena()->alloc(sizeof(PROB::ProbClosures));
  CNS::d_prob_parm =
      (PROB::ProbParm*)The_Arena()->alloc(sizeof(PROB::ProbParm));
  CNS::d_phys_bc = (BCRec*)The_Arena()->alloc(sizeof(BCRec));
#else
  CNS::d_prob_closures = h_prob_closures;
  CNS::d_prob_parm = h_prob_parm;
  CNS::d_phys_bc = h_phys_bc;
#endif

  // Read input parameters
  read_params();

  // Independent (solved) variables and their boundary condition types
  bool state_data_extrap = false;
  bool store_in_checkpoint = true;
  desc_lst.addDescriptor(State_Type, IndexType::TheCellType(),
                         StateDescriptor::Point, NGHOST, NCONS, &lincc_interp,
                         state_data_extrap, store_in_checkpoint);
  // https://github.com/AMReX-Codes/amrex/issues/396

  Vector<BCRec> bcs(NCONS);
  Vector<std::string> name(NCONS);

  // Physical boundary conditions ////////////////////////////////////////////
  int cnt = 0;
  set_scalar_bc(bcs[cnt], h_phys_bc);
  name[cnt] = "density";

  cnt++;
  set_x_vel_bc(bcs[cnt], h_phys_bc);
  name[cnt] = "xmom";

  cnt++;
  set_y_vel_bc(bcs[cnt], h_phys_bc);
  name[cnt] = "ymom";

#if (AMREX_SPACEDIM == 3)
  cnt++;
  set_z_vel_bc(bcs[cnt], h_phys_bc);
  name[cnt] = "zmom";
#endif
  cnt++;
  set_scalar_bc(bcs[cnt], h_phys_bc);
  name[cnt] = "energy";

  // PROB::ConsBCs
  // PROB::StateVarNames

  // Boundary conditions
  StateDescriptor::BndryFunc bndryfunc(cns_bcfill);
  StateDescriptor::setBndryFuncThreadSafety(true);
  bndryfunc.setRunOnGPU(true);
  // applies bndry func to all variables in desc_lst starting from from 0.
  desc_lst.setComponent(State_Type, 0, name, bcs, bndryfunc);
  ////////////////////////////////////////////////////////////////////////////

  // Define derived quantities ///////////////////////////////////////////////
  num_state_data_types = desc_lst.size();
  // Pressure
  derive_lst.add("pressure", IndexType::TheCellType(), 1, derpres,
                 the_same_box);
  derive_lst.addComponent("pressure", desc_lst, State_Type, h_prob_closures->URHO, NCONS);

  // Temperature
  derive_lst.add("temperature", IndexType::TheCellType(), 1, dertemp,
                 the_same_box);
  derive_lst.addComponent("temperature", desc_lst, State_Type, h_prob_closures->URHO, NCONS);

  // Velocities
  derive_lst.add("x_velocity", amrex::IndexType::TheCellType(), 1, dervel,
                 the_same_box);
  derive_lst.addComponent("x_velocity", desc_lst, State_Type, h_prob_closures->URHO, 1);
  derive_lst.addComponent("x_velocity", desc_lst, State_Type, h_prob_closures->UMX, 1);

  derive_lst.add("y_velocity", amrex::IndexType::TheCellType(), 1, dervel,
                 the_same_box);
  derive_lst.addComponent("y_velocity", desc_lst, State_Type, h_prob_closures->URHO, 1);
  derive_lst.addComponent("y_velocity", desc_lst, State_Type, h_prob_closures->UMY, 1);

#if (AMREX_SPACEDIM == 3)
  derive_lst.add("z_velocity", amrex::IndexType::TheCellType(), 1, dervel,
                 the_same_box);
  derive_lst.addComponent("z_velocity", desc_lst, State_Type, h_prob_closures->URHO, 1);
  derive_lst.addComponent("z_velocity", desc_lst, State_Type, h_prob_closures->UMZ, 1);
#endif
}