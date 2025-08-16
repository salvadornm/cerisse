#include "prob.H"

extern "C" {
void amrex_probinit(const int* /*init*/, const int* /*name*/, const int* /*namelen*/,
                    const amrex::Real* problo, const amrex::Real* probhi)
{
  // Parse params
  amrex::Real Pr = 0.7;
  {
    amrex::ParmParse pp("prob");
    pp.query("H", CNS::h_prob_parm->H);
    pp.query("Re_b", CNS::h_prob_parm->Re_b);
    pp.query("M_b", CNS::h_prob_parm->M_b);
    pp.query("Re_tau", CNS::h_prob_parm->Re_tau);
    pp.query("u_tau", CNS::h_prob_parm->u_tau);
    pp.query("Tw", CNS::h_prob_parm->T_w);
    pp.query("muw", CNS::h_prob_parm->mu_w);
    pp.query("Pr", Pr);
  }

  // Wall quantities
  CNS::h_prob_parm->rho_w = CNS::h_prob_parm->Re_tau * CNS::h_prob_parm->mu_w /
                            CNS::h_prob_parm->u_tau / CNS::h_prob_parm->H;

  // Bulk quantities
  auto eos = pele::physics::PhysicsType::eos();
  amrex::Real csw;
  eos.RTY2Cs(CNS::h_prob_parm->rho_w, CNS::h_prob_parm->T_w,
             CNS::h_prob_parm->massfrac.begin(), csw);
  CNS::h_prob_parm->u_b = CNS::h_prob_parm->M_b * csw;
  CNS::h_prob_parm->rho_b = CNS::h_prob_parm->Re_b * CNS::h_prob_parm->mu_w /
                            CNS::h_prob_parm->u_b / CNS::h_prob_parm->H;

  // Forcing per unit density
  CNS::h_prob_parm->f_x = CNS::h_prob_parm->rho_w * CNS::h_prob_parm->u_tau *
                          CNS::h_prob_parm->u_tau / CNS::h_prob_parm->H /
                          CNS::h_prob_parm->rho_b;

  Gpu::copyAsync(Gpu::hostToDevice, CNS::h_prob_parm, CNS::h_prob_parm + 1,
                 CNS::d_prob_parm);
  Gpu::streamSynchronize();

  // Use Sutherland law to approximate mu/mu_w = (T/T_w)^0.7
  auto& trans_parm = CNS::trans_parms.host_trans_parm();
  trans_parm.viscosity_mu_ref = CNS::h_prob_parm->mu_w;
  trans_parm.viscosity_T_ref = CNS::h_prob_parm->T_w;
  trans_parm.viscosity_S =
    CNS::h_prob_parm->T_w * 0.29; // use this to approximate mu ~ T^0.7 in T/Tw between 1 and 1.25
  trans_parm.Prandtl_number = Pr;
  trans_parm.const_bulk_viscosity = 0.0;
  trans_parm.const_diffusivity = 0.0;
  CNS::trans_parms.sync_to_device();

  // Print problem param
  amrex::Print() << "==================================="
                 << "\nrho_b = " << CNS::h_prob_parm->rho_b
                 << "\nu_b = " << CNS::h_prob_parm->u_b
                 << "\nmu_w = " << CNS::h_prob_parm->mu_w
                 << "\nc_w = " << csw
                 << "\nrho_w = " << CNS::h_prob_parm->rho_w
                 << "\nu_tau = " << CNS::h_prob_parm->u_tau
                 << "\nf_x = " << CNS::h_prob_parm->f_x
                 << "\n===================================\n";
}
}

void CNS::fill_ext_src(int i, int j, int k, amrex::Real time,
                       amrex::GeometryData const& geomdata,
                       amrex::Array4<const amrex::Real> const& state,
                       amrex::Array4<amrex::Real> const& ext_src, ProbParm const& pp)
{
  for (int nf = 0; nf <= NUM_FIELD; ++nf) {
    ext_src(i, j, k, nf * NVAR + UMX) += state(i, j, k, nf * NVAR + URHO) * pp.f_x;
    ext_src(i, j, k, nf * NVAR + UEDEN) += state(i, j, k, nf * NVAR + UMX) * pp.f_x;
  }
}

// Enable P(no I)D controller for bulk velocity (mass flow rate) by setting
// `USE_FULL_PROB_POST_TIMESTEP = TRUE` in GNUMakefile. 
// Adjust the P and D values below if necessary.
#if USE_FULL_PROB_POST_TIMESTEP
void CNS::full_prob_post_timestep(int /*iteration*/)
{
  // Sum bulk velcoity
  int finest_level = parent->finestLevel();
  amrex::Real bulk_u = 0.0;

  if (level == 0) {
    for (int lev = 0; lev <= finest_level; lev++) {
      CNS& cns_lev = getLevel(lev);

      amrex::MultiFab& S_new = cns_lev.get_new_data(State_Type);
      amrex::iMultiFab ifine_mask(cns_lev.grids, cns_lev.dmap, 1, 0);
      if (lev < parent->finestLevel()) {
        // mask out fine covered cells, do not sum
        ifine_mask =
          makeFineMask(cns_lev.grids, cns_lev.dmap, parent->boxArray(lev + 1),
                       cns_lev.fine_ratio, 1, 0);
      } else {
        ifine_mask.setVal(1);
      }
      amrex::MultiFab volume(cns_lev.grids, cns_lev.dmap, 1, 0);
      cns_lev.geom.GetVolume(volume);
      
      // const auto geomdata = cns_lev.geom.data();
      auto const& sarrs = S_new.const_arrays();
      auto const& marrs = ifine_mask.const_arrays();
      auto const& volarr = volume.const_arrays();

      auto reduce_tuple = amrex::ParReduce(
        TypeList<ReduceOpSum>{}, TypeList<Real>{}, S_new, IntVect(0),
        [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k) -> Real {
          const amrex::Real ux = sarrs[box_no](i, j, k, UMX) / sarrs[box_no](i, j, k, URHO);
          const amrex::Real mask = amrex::Real(marrs[box_no](i, j, k));
          const amrex::Real vol = volarr[box_no](i, j, k);
          return mask * vol * ux;
        });
      bulk_u += reduce_tuple;
    }
    // Reduction bulk_u
    amrex::ParallelDescriptor::ReduceRealSum(
      &bulk_u, 1, amrex::ParallelDescriptor::IOProcessorNumber());
    bulk_u /= (8.214 * 2.0 * 0.6845 * 2.0 * 1.369);  // / vol

    // Compute fx
    constexpr amrex::Real P = 0.1, D = 0.05;  // controller parameter
    const amrex::Real newtime = state[State_Type].curTime();
    const amrex::Real oldtime = state[State_Type].prevTime();
    amrex::Real dfx;
    if (amrex::ParallelDescriptor::IOProcessor()) {
      dfx = P * (CNS::d_prob_parm->u_b - bulk_u) / (newtime - oldtime) 
          + D * (CNS::d_prob_parm->bulk_u - bulk_u) / (newtime - oldtime);
    }
    amrex::ParallelDescriptor::Bcast(
      &dfx, 1, amrex::ParallelDescriptor::IOProcessorNumber());

    CNS::d_prob_parm->f_x += dfx;
    CNS::d_prob_parm->bulk_u = bulk_u;

    if (amrex::ParallelDescriptor::IOProcessor()) {
      // Write the quantities at this time
      const int log_index = 0;
      std::ostream& data_log = parent->DataLog(log_index);
      const int datwidth = 14;
      const int datprecision = 6;
      data_log << std::setw(datwidth) << newtime;
      data_log << std::setw(datwidth) << std::setprecision(datprecision) << bulk_u;
      data_log << std::setw(datwidth) << std::setprecision(datprecision) << CNS::d_prob_parm->f_x;
      data_log << std::endl;
    }
  }
}
#endif