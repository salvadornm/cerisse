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
//  (0)         (1)     (2)       (3)      (4)        (5)         (6)
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

static void set_z_vel_bc(BCRec& bc, const BCRec* phys_bc) {
  const int* lo_bc = phys_bc->lo();
  const int* hi_bc = phys_bc->hi();

  bc.setLo(0, tang_vel_bc[lo_bc[0]]);
  bc.setHi(0, tang_vel_bc[hi_bc[0]]);
#if (AMREX_SPACEDIM >= 2)  
  bc.setLo(1, tang_vel_bc[lo_bc[1]]);
  bc.setHi(1, tang_vel_bc[hi_bc[1]]);
#endif  
#if (AMREX_SPACEDIM == 3)  
  bc.setLo(2, norm_vel_bc[lo_bc[2]]);
  bc.setHi(2, norm_vel_bc[hi_bc[2]]);
#endif  

}

void CNS::variableSetUp() {

  amrex::Print( ) << " oo CNS::variableSetUp " << std::endl; 

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
                         StateDescriptor::Point, h_prob_closures->NGHOST, h_prob_closures->NCONS, &lincc_interp,
                         state_data_extrap, store_in_checkpoint);
  // https://github.com/AMReX-Codes/amrex/issues/396

  // Setup Stats_Type
  if (h_prob_closures->NSTAT > 0 ) {
    compute_stats=true;
    desc_lst.addDescriptor(Stats_Type, IndexType::TheCellType(),
                          StateDescriptor::Point, 0, h_prob_closures->NSTAT, &lincc_interp,
                          state_data_extrap, store_in_checkpoint);
  }                        
  
  // Physical boundary conditions ////////////////////////////////////////////
  Vector<BCRec> bcs(PROB::ProbClosures::NCONS);
  Vector<int> cons_vars_type = indicies_t::get_cons_vars_type();

  for (int cnt=0;cnt<h_prob_closures->NCONS;cnt++) {

    switch (cons_vars_type[cnt])
    {
      case 0:
        set_scalar_bc(bcs[cnt], h_phys_bc);
        break;
      case 1:
        set_x_vel_bc (bcs[cnt], h_phys_bc);
        break;
      case 2:  
        set_y_vel_bc (bcs[cnt], h_phys_bc); 
        break;
      case 3:  
        set_z_vel_bc (bcs[cnt], h_phys_bc);    
        break;
      default:
        std::cout << " error ... variableSetUp" << std::endl;
        exit(1); 
    }

  }

  // Boundary conditions
  StateDescriptor::BndryFunc bndryfunc(cns_bcfill);
  StateDescriptor::setBndryFuncThreadSafety(true);
  bndryfunc.setRunOnGPU(true);
  // applies bndry func to all variables in desc_lst starting from from 0.
  desc_lst.setComponent(State_Type, 0, PROB::ProbClosures::get_cons_vars_names(), bcs, bndryfunc);

  num_state_data_types = desc_lst.size();

  //printf("num_state_data_types (1) %d \n",num_state_data_types);
  // exit(0);

  // SET-UP Stats Type 
  ////////////////////////////////////////////////////////////////////////////
  if (h_prob_closures->NSTAT > 0) { 

    int NSTAT     = h_prob_closures->NSTAT; 
    int NSTAT_VEL = h_prob_closures->NSTAT_VEL;
    
    Vector<BCRec> stats_bcs(NSTAT);
    Vector<std::string>  stats_name(NSTAT);
  
    // names velocity
    if (NSTAT_VEL > 0) {
      int statv = 0;
      stats_name[statv] = "x_velocityMEAN";
#if AMREX_SPACEDIM >1    
      statv++; stats_name[statv] = "y_velocityMEAN";
#endif
#if AMREX_SPACEDIM == 3    
      statv++; stats_name[statv] = "z_velocityMEAN";
#endif
      statv++; stats_name[statv] = "x_velocitySQR";
#if AMREX_SPACEDIM >1    
      statv++; stats_name[statv] = "y_velocitySQR";
#endif
#if AMREX_SPACEDIM == 3    
      statv++; stats_name[statv] = "z_velocitySQR";
#endif
#if AMREX_SPACEDIM >1    
      statv++; stats_name[statv] = "xy_velocityMEAN";
#endif
#if AMREX_SPACEDIM == 3    
      statv++; stats_name[statv] = "xz_velocityMEAN";
      statv++; stats_name[statv] = "yz_velocityMEAN";    
#endif  
    
    }

    // bc
    for (int statv=0;statv<h_prob_closures->NSTAT;statv++) {
      set_scalar_bc(stats_bcs[statv], h_phys_bc);
    }
    StateDescriptor::BndryFunc bndryfuncstats( cns_bcfill);
    bndryfuncstats.setRunOnGPU(true);
    desc_lst.setComponent(Stats_Type, 0, stats_name, stats_bcs, bndryfuncstats);
  }
  //printf("num_state_data_types(2) %d \n",num_state_data_types);

  ////////////////////////////////////////////////////////////////////////////

  // Define derived quantities ///////////////////////////////////////////////
  // Pressure
  derive_lst.add("pressure", IndexType::TheCellType(), 1, derpres,
                 the_same_box);
  derive_lst.addComponent("pressure", desc_lst, State_Type, 0, h_prob_closures->NCONS);

  // Temperature
  derive_lst.add("temperature", IndexType::TheCellType(), 1, dertemp,
                 the_same_box);
  derive_lst.addComponent("temperature", desc_lst, State_Type, 0, h_prob_closures->NCONS);

  // Velocities
  derive_lst.add("velocity", IndexType::TheCellType(), amrex::SpaceDim,
                 {AMREX_D_DECL("x_velocity", "y_velocity", "z_velocity")},
                 dervel, the_same_box);
  derive_lst.addComponent("velocity", desc_lst, State_Type, 0,
                          h_prob_closures->NCONS);

  // Density
  derive_lst.add("density", IndexType::TheCellType(), 1, derdensity,
                 the_same_box);
  derive_lst.addComponent("density", desc_lst, State_Type, 0, h_prob_closures->NCONS);

  // Kinetic energy
  derive_lst.add("kinetic_energy", IndexType::TheCellType(), 1, derkineticenergy,
                 the_same_box);
  derive_lst.addComponent("kinetic_energy", desc_lst, State_Type, 0, h_prob_closures->NCONS);

  // Vorticity magnitude
  derive_lst.add("magvort", IndexType::TheCellType(), 1, dermagvort,
                 DeriveRec::GrowBoxByOne);
  derive_lst.addComponent("magvort", desc_lst, State_Type, 0, h_prob_closures->NCONS);

  // Enstrophy
  derive_lst.add("enstrophy", IndexType::TheCellType(), 1, derenstrophy,
                 DeriveRec::GrowBoxByOne);
  derive_lst.addComponent("enstrophy", desc_lst, State_Type, 0, h_prob_closures->NCONS);
}