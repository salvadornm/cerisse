#include "CNS.H"
#include "derive.H"
#include "prob.H"

using namespace amrex;

int CNS::num_state_data_types = 0;
Parm* CNS::h_parm = nullptr;
Parm* CNS::d_parm = nullptr;
ProbParm* CNS::h_prob_parm = nullptr;
ProbParm* CNS::d_prob_parm = nullptr;

static Box the_same_box (const Box& b) { return b; }
//static Box grow_box_by_one (const Box& b) { return amrex::grow(b,1); }

using BndryFunc = StateDescriptor::BndryFunc;

// Components are:
//  Interior,        Inflow,          Outflow,          Symmetry,             SlipWall,             NoSlipWall,           UserBC
static int scalar_bc[] = {
    BCType::int_dir, BCType::ext_dir, BCType::foextrap, BCType::reflect_even, BCType::reflect_even, BCType::reflect_even, BCType::ext_dir
};

static int norm_vel_bc[] = {
    BCType::int_dir, BCType::ext_dir, BCType::foextrap, BCType::reflect_odd,  BCType::reflect_odd,  BCType::reflect_odd,  BCType::ext_dir
};

static int tang_vel_bc[] = {
    BCType::int_dir, BCType::ext_dir, BCType::foextrap, BCType::reflect_even, BCType::reflect_even, BCType::reflect_odd,  BCType::ext_dir
};

static int react_src_bc[] = {
    BCType::int_dir, BCType::reflect_even, BCType::reflect_even, BCType::reflect_even, BCType::reflect_even, BCType::reflect_even, BCType::reflect_even
};

static
void
set_scalar_bc(BCRec& bc, const BCRec& phys_bc)
{
    const int* lo_bc = phys_bc.lo();
    const int* hi_bc = phys_bc.hi();
    for (int dir = 0; dir < AMREX_SPACEDIM; dir++)
    {
        bc.setLo(dir, scalar_bc[lo_bc[dir]]);
        bc.setHi(dir, scalar_bc[hi_bc[dir]]);
    }
}

static
void
set_x_vel_bc(BCRec& bc, const BCRec& phys_bc)
{
    const int* lo_bc = phys_bc.lo();
    const int* hi_bc = phys_bc.hi();
    AMREX_D_TERM(
        bc.setLo(0, norm_vel_bc[lo_bc[0]]);  bc.setHi(0, norm_vel_bc[hi_bc[0]]); ,
        bc.setLo(1, tang_vel_bc[lo_bc[1]]);  bc.setHi(1, tang_vel_bc[hi_bc[1]]); ,
        bc.setLo(2, tang_vel_bc[lo_bc[2]]);  bc.setHi(2, tang_vel_bc[hi_bc[2]]); 
    )
}

static
void
set_y_vel_bc(BCRec& bc, const BCRec& phys_bc)
{
    const int* lo_bc = phys_bc.lo();
    const int* hi_bc = phys_bc.hi();
    AMREX_D_TERM(
        bc.setLo(0, tang_vel_bc[lo_bc[0]]);  bc.setHi(0, tang_vel_bc[hi_bc[0]]); ,
        bc.setLo(1, norm_vel_bc[lo_bc[1]]);  bc.setHi(1, norm_vel_bc[hi_bc[1]]); ,
        bc.setLo(2, tang_vel_bc[lo_bc[2]]);  bc.setHi(2, tang_vel_bc[hi_bc[2]]); 
    )
}

#if (AMREX_SPACEDIM == 3)
static
void
set_z_vel_bc(BCRec& bc, const BCRec& phys_bc)
{
    const int* lo_bc = phys_bc.lo();
    const int* hi_bc = phys_bc.hi();
    bc.setLo(0,tang_vel_bc[lo_bc[0]]);  bc.setHi(0,tang_vel_bc[hi_bc[0]]);
    bc.setLo(1,tang_vel_bc[lo_bc[1]]);  bc.setHi(1,tang_vel_bc[hi_bc[1]]);
    bc.setLo(2,norm_vel_bc[lo_bc[2]]);  bc.setHi(2,norm_vel_bc[hi_bc[2]]);
}
#endif

static 
void
set_react_src_bc(BCRec& bc, const BCRec& phys_bc)
{
    const int* lo_bc = phys_bc.lo();
    const int* hi_bc = phys_bc.hi();
    for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
        bc.setLo(dir, react_src_bc[lo_bc[dir]]);
        bc.setHi(dir, react_src_bc[hi_bc[dir]]);
    }
}

void
CNS::variableSetUp()
{
    h_parm = new Parm{}; // This is deleted in CNS::variableCleanUp().
    h_prob_parm = new ProbParm{};
    d_parm = (Parm*)The_Arena()->alloc(sizeof(Parm)); //static_cast
    d_prob_parm = (ProbParm*)The_Arena()->alloc(sizeof(ProbParm));
    trans_parms.allocate();

    read_params();

    // Setup State_Type
    bool state_data_extrap = false;
    bool store_in_checkpoint = true;
    desc_lst.addDescriptor(State_Type, IndexType::TheCellType(),
                           StateDescriptor::Point, NUM_GROW, NUM_STATE,
                           &eb_mf_cell_cons_interp, state_data_extrap, store_in_checkpoint);

    amrex::Vector<amrex::BCRec> bcs(NUM_STATE);
    amrex::Vector<std::string>  name(NUM_STATE);
    amrex::BCRec bc;
    int cnt = 0;
           set_scalar_bc(bc,phys_bc); bcs[cnt] = bc; name[cnt] = "density";
    cnt++; set_x_vel_bc(bc,phys_bc);  bcs[cnt] = bc; name[cnt] = "xmom";
    cnt++; set_y_vel_bc(bc,phys_bc);  bcs[cnt] = bc; name[cnt] = "ymom";
#if (AMREX_SPACEDIM == 3)
    cnt++; set_z_vel_bc(bc,phys_bc);  bcs[cnt] = bc; name[cnt] = "zmom";
#endif
    cnt++; set_scalar_bc(bc,phys_bc); bcs[cnt] = bc; name[cnt] = "rho_E";
    cnt++; set_scalar_bc(bc,phys_bc); bcs[cnt] = bc; name[cnt] = "rho_e";
    cnt++; set_scalar_bc(bc,phys_bc); bcs[cnt] = bc; name[cnt] = "Temp";
    
    pele::physics::eos::speciesNames<
        pele::physics::PhysicsType::eos_type>(spec_names);

    if (amrex::ParallelDescriptor::IOProcessor()) {
        Print() << NUM_SPECIES << " Species: ";
        for (int i = 0; i < NUM_SPECIES; i++) {
            Print() << spec_names[i] << ' ' << ' ';
        }
        Print() << std::endl;
    }
    for (int i = 0; i < NUM_SPECIES; ++i) {
        cnt++; 
        set_scalar_bc(bc, phys_bc); 
        bcs[cnt] = bc; 
        name[cnt] = "rho_" + spec_names[i];
    }

    StateDescriptor::BndryFunc bndryfunc(cns_bcfill);
    bndryfunc.setRunOnGPU(true);

    desc_lst.setComponent(State_Type, URHO, name, bcs, bndryfunc);

    // Setup React_Type
    store_in_checkpoint = do_react;
    desc_lst.addDescriptor(Reactions_Type, amrex::IndexType::TheCellType(),
                           amrex::StateDescriptor::Point, 0, NUM_SPECIES+2, 
                           &eb_mf_cell_cons_interp, state_data_extrap, store_in_checkpoint);

    amrex::Vector<amrex::BCRec> react_bcs(NUM_SPECIES+2);
    amrex::Vector<std::string>  react_name(NUM_SPECIES+2);

    for (int i = 0; i < NUM_SPECIES; ++i) {
        set_react_src_bc(bc, phys_bc); 
        react_bcs[i] = bc; 
        react_name[i] = "rho_omega_" + spec_names[i];
    }
    set_react_src_bc(bc, phys_bc);
    react_bcs[NUM_SPECIES] = bc;
    react_name[NUM_SPECIES] = "rhoe_dot";
    set_react_src_bc(bc, phys_bc);
    react_bcs[NUM_SPECIES + 1] = bc;
    react_name[NUM_SPECIES + 1] = "heatRelease";

    amrex::StateDescriptor::BndryFunc bndryfunc2(cns_react_bcfill);
    bndryfunc2.setRunOnGPU(true);

    desc_lst.setComponent(Reactions_Type, 0, react_name, react_bcs, bndryfunc2);

    // Setup Cost_Type
    desc_lst.addDescriptor(Cost_Type, IndexType::TheCellType(), StateDescriptor::Point,
                           0, 1, &pc_interp);
    desc_lst.setComponent(Cost_Type, 0, "Cost", bc, bndryfunc);

    num_state_data_types = desc_lst.size();

    StateDescriptor::setBndryFuncThreadSafety(true);

    // DEFINE DERIVED QUANTITIES
    // Pressure
    derive_lst.add("pressure",IndexType::TheCellType(),1,
                   cns_derpres,the_same_box);
    derive_lst.addComponent("pressure",desc_lst,State_Type,UEINT,1);

    // Velocities
    derive_lst.add("x_velocity",IndexType::TheCellType(),1,
                   cns_dervel,the_same_box);
    derive_lst.addComponent("x_velocity",desc_lst,State_Type,URHO,1);
    derive_lst.addComponent("x_velocity",desc_lst,State_Type,UMX,1);

    derive_lst.add("y_velocity",IndexType::TheCellType(),1,
                   cns_dervel,the_same_box);
    derive_lst.addComponent("y_velocity",desc_lst,State_Type,URHO,1);
    derive_lst.addComponent("y_velocity",desc_lst,State_Type,UMY,1);

#if (AMREX_SPACEDIM == 3)
    derive_lst.add("z_velocity",IndexType::TheCellType(),1,
                   cns_dervel,the_same_box);
    derive_lst.addComponent("z_velocity",desc_lst,State_Type,URHO,1);
    derive_lst.addComponent("z_velocity",desc_lst,State_Type,UMZ,1);
#endif
}

void
CNS::variableCleanUp()
{
    delete h_parm;
    delete h_prob_parm;
    The_Arena()->free(d_parm);
    The_Arena()->free(d_prob_parm);
    desc_lst.clear();
    derive_lst.clear();

#ifdef AMREX_USE_GPU
    The_Arena()->free(dp_refine_boxes);
#endif
}
