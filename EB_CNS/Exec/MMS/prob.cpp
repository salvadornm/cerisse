#include "prob.H"

extern "C" {
void amrex_probinit(const int* /*init*/, const int* /*name*/, const int* /*namelen*/,
                    const amrex::Real* problo, const amrex::Real* probhi)
{
  // Parse params
  bool do_visc = false;
  {
    amrex::ParmParse pp("prob");
    pp.query("reynolds", CNS::h_prob_parm->reynolds);
    pp.query("mach", CNS::h_prob_parm->mach);
    pp.query("prandtl", CNS::h_prob_parm->prandtl);
    pp.query("p0", CNS::h_prob_parm->p0);
    pp.query("T0", CNS::h_prob_parm->T0);
    pp.query("rho0", CNS::h_prob_parm->rho0);
    pp.query("rho_x_fact", CNS::h_prob_parm->rho_x_fact);
    pp.query("rho_y_fact", CNS::h_prob_parm->rho_y_fact);
    pp.query("rho_z_fact", CNS::h_prob_parm->rho_z_fact);
    pp.query("u_x_fact", CNS::h_prob_parm->u_x_fact);
    pp.query("u_y_fact", CNS::h_prob_parm->u_y_fact);
    pp.query("u_z_fact", CNS::h_prob_parm->u_z_fact);
    pp.query("v_0_fact", CNS::h_prob_parm->v_0_fact);
    pp.query("v_x_fact", CNS::h_prob_parm->v_x_fact);
    pp.query("v_y_fact", CNS::h_prob_parm->v_y_fact);
    pp.query("v_z_fact", CNS::h_prob_parm->v_z_fact);
    pp.query("w_0_fact", CNS::h_prob_parm->w_0_fact);
    pp.query("w_x_fact", CNS::h_prob_parm->w_x_fact);
    pp.query("w_y_fact", CNS::h_prob_parm->w_y_fact);
    pp.query("w_z_fact", CNS::h_prob_parm->w_z_fact);
    pp.query("p_x_fact", CNS::h_prob_parm->p_x_fact);
    pp.query("p_y_fact", CNS::h_prob_parm->p_y_fact);
    pp.query("p_z_fact", CNS::h_prob_parm->p_z_fact);
    pp.query("a_rhox", CNS::h_prob_parm->a_rhox);
    pp.query("a_rhoy", CNS::h_prob_parm->a_rhoy);
    pp.query("a_rhoz", CNS::h_prob_parm->a_rhoz);
    pp.query("a_ux", CNS::h_prob_parm->a_ux);
    pp.query("a_uy", CNS::h_prob_parm->a_uy);
    pp.query("a_uz", CNS::h_prob_parm->a_uz);
    pp.query("a_vx", CNS::h_prob_parm->a_vx);
    pp.query("a_vy", CNS::h_prob_parm->a_vy);
    pp.query("a_vz", CNS::h_prob_parm->a_vz);
    pp.query("a_wx", CNS::h_prob_parm->a_wx);
    pp.query("a_wy", CNS::h_prob_parm->a_wy);
    pp.query("a_wz", CNS::h_prob_parm->a_wz);
    pp.query("a_px", CNS::h_prob_parm->a_px);
    pp.query("a_py", CNS::h_prob_parm->a_py);
    pp.query("a_pz", CNS::h_prob_parm->a_pz);
    // pp.query("masa_solution_name", masa_solution_name);

    amrex::ParmParse ppc("cns");
    ppc.query("do_visc", do_visc);
  }

  // Define the length scale
  CNS::h_prob_parm->L = probhi[0] - problo[0]; // L_x

  // Initial density, velocity, and material properties
  amrex::Real eint;
  amrex::Real cs;
  amrex::Real cp;
  amrex::Real massfrac[NUM_SPECIES] = {1.0};
  auto eos = pele::physics::PhysicsType::eos();
  eos.PYT2RE(CNS::h_prob_parm->p0, massfrac, CNS::h_prob_parm->T0,
             CNS::h_prob_parm->rho0, eint);
  eos.RTY2Cs(CNS::h_prob_parm->rho0, CNS::h_prob_parm->T0, massfrac, cs);
  eos.TY2Cp(CNS::h_prob_parm->T0, massfrac, cp);
  CNS::h_prob_parm->u0 = CNS::h_prob_parm->mach * cs;

  auto& trans_parm = CNS::trans_parms.host_trans_parm();
  if (do_visc) {
    trans_parm.const_bulk_viscosity = 0.0;
    trans_parm.const_diffusivity = 0.0;
    trans_parm.const_viscosity = CNS::h_prob_parm->rho0 * CNS::h_prob_parm->u0 *
                                 CNS::h_prob_parm->L / CNS::h_prob_parm->reynolds;
    trans_parm.const_conductivity =
      trans_parm.const_viscosity * cp / CNS::h_prob_parm->prandtl;
  } else {
    trans_parm.const_bulk_viscosity = 0.0;
    trans_parm.const_diffusivity = 0.0;
    trans_parm.const_viscosity = 0.0;
    trans_parm.const_conductivity = 0.0;
  }
  CNS::trans_parms.sync_to_device();

  Gpu::copy(Gpu::hostToDevice, CNS::h_prob_parm, CNS::h_prob_parm + 1,
            CNS::d_prob_parm);

  // MASA parameters for the following functions (masa/src/cns.cpp)
  // rho = rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L)
  // u = u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L)
  // ...
  // p = p_0 + p_x * cos(a_px * pi * x / L) + p_y * sin(a_py * pi * y / L) + p_z * cos(a_pz * pi * z / L)
  masa_init("mms", "navierstokes_3d_compressible");
  masa_set_param("L", CNS::h_prob_parm->L);
  masa_set_param("R",
                 pele::physics::Constants::RU / pele::physics::Constants::AIRMW);
  masa_set_param("k", trans_parm.const_conductivity);
  masa_set_param("Gamma", eos.gamma);
  masa_set_param("mu", trans_parm.const_viscosity);
  masa_set_param("rho_0", CNS::h_prob_parm->rho0);
  masa_set_param("rho_x", CNS::h_prob_parm->rho_x_fact * CNS::h_prob_parm->rho0);
  masa_set_param("rho_y", CNS::h_prob_parm->rho_y_fact * CNS::h_prob_parm->rho0);
  masa_set_param("rho_z", CNS::h_prob_parm->rho_z_fact * CNS::h_prob_parm->rho0);
  masa_set_param("u_0", CNS::h_prob_parm->u0);
  masa_set_param("u_x", CNS::h_prob_parm->u_x_fact * CNS::h_prob_parm->u0);
  masa_set_param("u_y", CNS::h_prob_parm->u_y_fact * CNS::h_prob_parm->u0);
  masa_set_param("u_z", CNS::h_prob_parm->u_z_fact * CNS::h_prob_parm->u0);
  masa_set_param("v_0", CNS::h_prob_parm->v_0_fact * CNS::h_prob_parm->u0);
  masa_set_param("v_x", CNS::h_prob_parm->v_x_fact * CNS::h_prob_parm->u0);
  masa_set_param("v_y", CNS::h_prob_parm->v_y_fact * CNS::h_prob_parm->u0);
  masa_set_param("v_z", CNS::h_prob_parm->v_z_fact * CNS::h_prob_parm->u0);
  masa_set_param("w_0", CNS::h_prob_parm->w_0_fact * CNS::h_prob_parm->u0);
  masa_set_param("w_x", CNS::h_prob_parm->w_x_fact * CNS::h_prob_parm->u0);
  masa_set_param("w_y", CNS::h_prob_parm->w_y_fact * CNS::h_prob_parm->u0);
  masa_set_param("w_z", CNS::h_prob_parm->w_z_fact * CNS::h_prob_parm->u0);
  masa_set_param("p_0", CNS::h_prob_parm->p0);
  masa_set_param("p_x", CNS::h_prob_parm->p_x_fact * CNS::h_prob_parm->p0);
  masa_set_param("p_y", CNS::h_prob_parm->p_y_fact * CNS::h_prob_parm->p0);
  masa_set_param("p_z", CNS::h_prob_parm->p_z_fact * CNS::h_prob_parm->p0);
  masa_set_param("a_rhox", CNS::h_prob_parm->a_rhox);
  masa_set_param("a_rhoy", CNS::h_prob_parm->a_rhoy);
  masa_set_param("a_rhoz", CNS::h_prob_parm->a_rhoz);
  masa_set_param("a_ux", CNS::h_prob_parm->a_ux);
  masa_set_param("a_uy", CNS::h_prob_parm->a_uy);
  masa_set_param("a_uz", CNS::h_prob_parm->a_uz);
  masa_set_param("a_vx", CNS::h_prob_parm->a_vx);
  masa_set_param("a_vy", CNS::h_prob_parm->a_vy);
  masa_set_param("a_vz", CNS::h_prob_parm->a_vz);
  masa_set_param("a_wx", CNS::h_prob_parm->a_wx);
  masa_set_param("a_wy", CNS::h_prob_parm->a_wy);
  masa_set_param("a_wz", CNS::h_prob_parm->a_wz);
  masa_set_param("a_px", CNS::h_prob_parm->a_px);
  masa_set_param("a_py", CNS::h_prob_parm->a_py);
  masa_set_param("a_pz", CNS::h_prob_parm->a_pz);

  if (amrex::ParallelDescriptor::IOProcessor()) { masa_display_param(); }
  masa_sanity_check();
}
}

void CNS::fill_ext_src(int i, int j, int k, amrex::Real time,
                       amrex::GeometryData const& geomdata,
                       amrex::Array4<const amrex::Real> const& state,
                       amrex::Array4<amrex::Real> const& ext_src,
                       ProbParm const& /*prob_parm*/)
{
  // Geometry
  const amrex::Real* prob_lo = geomdata.ProbLo();
  const amrex::Real* dx = geomdata.CellSize();
  const amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];
  const amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1];
  const amrex::Real z = prob_lo[2] + (k + 0.5) * dx[2];

  // Get source term from MASA
  amrex::Real source_term[5] = {0.0};  
  source_term[0] = masa_eval_3d_source_rho(x, y, z);
  source_term[1] = masa_eval_3d_source_rho_u(x, y, z);
  source_term[2] = masa_eval_3d_source_rho_v(x, y, z);
  source_term[3] = masa_eval_3d_source_rho_w(x, y, z);
  source_term[4] = masa_eval_3d_source_rho_e(x, y, z);

  // Set the state
  for (int nf = 0; nf <= NUM_FIELD; ++nf) {
    ext_src(i, j, k, nf * NVAR + URHO) += source_term[0];
    ext_src(i, j, k, nf * NVAR + UMX) += source_term[1];
    ext_src(i, j, k, nf * NVAR + UMY) += source_term[2];
    ext_src(i, j, k, nf * NVAR + UMZ) += source_term[3];
    ext_src(i, j, k, nf * NVAR + UEDEN) += source_term[4];
    amrex::Real rhoinv = 1.0 / state(i, j, k, nf * NVAR + URHO);
    amrex::Real Y;
    for (int n = 0; n < NUM_SPECIES; n++) {
      Y = state(i, j, k, nf * NVAR + UFS + n) * rhoinv;
      ext_src(i, j, k, nf * NVAR + UFS + n) += source_term[0] * Y;
    }
  }
}

void CNS::full_prob_post_timestep(int /*iteration*/)
{
  int finest_level = parent->finestLevel();
  amrex::Real time = state[State_Type].curTime();
  amrex::Real rho_mms_err = 0.0;
  amrex::Real u_mms_err = 0.0;
  amrex::Real v_mms_err = 0.0;
  amrex::Real w_mms_err = 0.0;
  amrex::Real p_mms_err = 0.0;
  amrex::Real rho_residual = 0.0;
  amrex::Real rhou_residual = 0.0;
  amrex::Real rhov_residual = 0.0;
  amrex::Real rhow_residual = 0.0;
  amrex::Real rhoE_residual = 0.0;

  if (level == 0) {
    if (amrex::ParallelDescriptor::IOProcessor()) {
      amrex::Print() << "... MMS problem post timestep" << std::endl;
    }

    // Calculate the errors and residuals
    for (int lev = 0; lev <= finest_level; lev++) {
      CNS& cns_lev = getLevel(lev);

      amrex::MultiFab& S_new = cns_lev.get_new_data(State_Type);
      amrex::MultiFab& S_old = cns_lev.get_old_data(State_Type);
      amrex::iMultiFab ifine_mask(cns_lev.grids, cns_lev.dmap, 1, 0);
      if (lev < parent->finestLevel()) {
        ifine_mask =
          makeFineMask(cns_lev.grids, cns_lev.dmap, parent->boxArray(lev + 1),
                       cns_lev.fine_ratio, 1, 0);
      } else {
        ifine_mask.setVal(1);
      }

      const auto geomdata = cns_lev.geom.data();
      auto const& sarrs = S_new.const_arrays();
      auto const& soldarrs = S_old.const_arrays();
      auto const& marrs = ifine_mask.const_arrays();

      auto reduce_tuple = amrex::ParReduce(
        TypeList<ReduceOpSum, ReduceOpSum, ReduceOpSum, ReduceOpSum, ReduceOpSum,
                 ReduceOpSum, ReduceOpSum, ReduceOpSum, ReduceOpSum, ReduceOpSum>{},
        TypeList<Real, Real, Real, Real, Real, Real, Real, Real, Real, Real>{},
        S_new, IntVect(0),
        [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k)
          -> GpuTuple<Real, Real, Real, Real, Real, Real, Real, Real, Real, Real> {
          const amrex::Real* prob_lo = geomdata.ProbLo();
          const amrex::Real* dx = geomdata.CellSize();
          const amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];
          const amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1];
          const amrex::Real z = prob_lo[2] + (k + 0.5) * dx[2];

          // Get exact solution (point value)
          const amrex::Real rho_exact = masa_eval_3d_exact_rho(x, y, z);
          const amrex::Real u_exact = masa_eval_3d_exact_u(x, y, z);
          const amrex::Real v_exact = masa_eval_3d_exact_v(x, y, z);
          const amrex::Real w_exact = masa_eval_3d_exact_w(x, y, z);
          const amrex::Real p_exact = masa_eval_3d_exact_p(x, y, z);
                    
          // Get numerical solution
          const amrex::Real rho = sarrs[box_no](i, j, k, URHO);
          const amrex::Real rhoinv = 1.0 / rho;
          const amrex::Real u = sarrs[box_no](i, j, k, UMX) * rhoinv;
          const amrex::Real v = sarrs[box_no](i, j, k, UMY) * rhoinv;
          const amrex::Real w = sarrs[box_no](i, j, k, UMZ) * rhoinv;
          const amrex::Real ei =
            sarrs[box_no](i, j, k, UEDEN) * rhoinv - 0.5 * (u * u + v * v + w * w);
          amrex::Real Y[NUM_SPECIES];
          for (int n = 0; n < NUM_SPECIES; n++)
            Y[n] = sarrs[box_no](i, j, k, UFS + n) * rhoinv;
          amrex::Real T, p;
          auto eos = pele::physics::PhysicsType::eos();
          eos.REY2T(rho, ei, Y, T);
          eos.RTY2P(rho, T, Y, p);

          // Calculate errors
          const amrex::Real mask =
            amrex::Real(marrs[box_no](i, j, k)); // mask out fine-covered cells
          const amrex::Real vol = amrex::Geometry::Volume({i, j, k}, geomdata);
          amrex::Real rho_err = mask * vol * (rho - rho_exact) * (rho - rho_exact);
          amrex::Real u_err = mask * vol * (u - u_exact) * (u - u_exact);
          amrex::Real v_err = mask * vol * (v - v_exact) * (v - v_exact);
          amrex::Real w_err = mask * vol * (w - w_exact) * (w - w_exact);
          amrex::Real p_err = mask * vol * (p - p_exact) * (p - p_exact);

          // Calculate residuals
          const amrex::Real ru = sarrs[box_no](i, j, k, UMX);
          const amrex::Real rv = sarrs[box_no](i, j, k, UMY);
          const amrex::Real rw = sarrs[box_no](i, j, k, UMZ);
          const amrex::Real re = sarrs[box_no](i, j, k, UEDEN);
          const amrex::Real rho_old = soldarrs[box_no](i, j, k, URHO);
          const amrex::Real ru_old = soldarrs[box_no](i, j, k, UMX);
          const amrex::Real rv_old = soldarrs[box_no](i, j, k, UMY);
          const amrex::Real rw_old = soldarrs[box_no](i, j, k, UMZ);
          const amrex::Real re_old = soldarrs[box_no](i, j, k, UEDEN);
          amrex::Real rho_res = mask * vol * (rho - rho_old) * (rho - rho_old);
          amrex::Real ru_res = mask * vol * (ru - ru_old) * (ru - ru_old);
          amrex::Real rv_res = mask * vol * (rv - rv_old) * (rv - rv_old);
          amrex::Real rw_res = mask * vol * (rw - rw_old) * (rw - rw_old);
          amrex::Real re_res = mask * vol * (re - re_old) * (re - re_old);

          return {rho_err, u_err,  v_err,  w_err,  p_err,
                  rho_res, ru_res, rv_res, rw_res, re_res};
        });

      rho_mms_err += amrex::get<0>(reduce_tuple);
      u_mms_err += amrex::get<1>(reduce_tuple);
      v_mms_err += amrex::get<2>(reduce_tuple);
      w_mms_err += amrex::get<3>(reduce_tuple);
      p_mms_err += amrex::get<4>(reduce_tuple);
      rho_residual += amrex::get<5>(reduce_tuple);
      rhou_residual += amrex::get<6>(reduce_tuple);
      rhov_residual += amrex::get<7>(reduce_tuple);
      rhow_residual += amrex::get<8>(reduce_tuple);
      rhoE_residual += amrex::get<9>(reduce_tuple);
    }

    // Reductions
    amrex::ParallelDescriptor::ReduceRealSum(
      &rho_mms_err, 1, amrex::ParallelDescriptor::IOProcessorNumber());
    amrex::ParallelDescriptor::ReduceRealSum(
      &u_mms_err, 1, amrex::ParallelDescriptor::IOProcessorNumber());
    amrex::ParallelDescriptor::ReduceRealSum(
      &v_mms_err, 1, amrex::ParallelDescriptor::IOProcessorNumber());
    amrex::ParallelDescriptor::ReduceRealSum(
      &w_mms_err, 1, amrex::ParallelDescriptor::IOProcessorNumber());
    amrex::ParallelDescriptor::ReduceRealSum(
      &p_mms_err, 1, amrex::ParallelDescriptor::IOProcessorNumber());
    amrex::ParallelDescriptor::ReduceRealSum(
      &rho_residual, 1, amrex::ParallelDescriptor::IOProcessorNumber());
    amrex::ParallelDescriptor::ReduceRealSum(
      &rhou_residual, 1, amrex::ParallelDescriptor::IOProcessorNumber());
    amrex::ParallelDescriptor::ReduceRealSum(
      &rhov_residual, 1, amrex::ParallelDescriptor::IOProcessorNumber());
    amrex::ParallelDescriptor::ReduceRealSum(
      &rhow_residual, 1, amrex::ParallelDescriptor::IOProcessorNumber());
    amrex::ParallelDescriptor::ReduceRealSum(
      &rhoE_residual, 1, amrex::ParallelDescriptor::IOProcessorNumber());

    // Get the norm
    rho_mms_err = std::sqrt(rho_mms_err);
    u_mms_err = std::sqrt(u_mms_err);
    v_mms_err = std::sqrt(v_mms_err);
    w_mms_err = std::sqrt(w_mms_err);
    p_mms_err = std::sqrt(p_mms_err);
    rho_residual = std::sqrt(rho_residual);
    rhou_residual = std::sqrt(rhou_residual);
    rhov_residual = std::sqrt(rhov_residual);
    rhow_residual = std::sqrt(rhow_residual);
    rhoE_residual = std::sqrt(rhoE_residual);

    if (amrex::ParallelDescriptor::IOProcessor()) {
      amrex::Print() << "TIME = " << time << '\n';
      amrex::Print() << " RHO MMS ERROR  = " << rho_mms_err << '\n';
      amrex::Print() << " U MMS ERROR    = " << u_mms_err << '\n';
      amrex::Print() << " V MMS ERROR    = " << v_mms_err << '\n';
      amrex::Print() << " W MMS ERROR    = " << w_mms_err << '\n';
      amrex::Print() << " P MMS ERROR    = " << p_mms_err << '\n';
      amrex::Print() << " RHO RESIDUAL   = " << rho_residual << '\n';
      amrex::Print() << " RHO*U RESIDUAL = " << rhou_residual << '\n';
      amrex::Print() << " RHO*V RESIDUAL = " << rhov_residual << '\n';
      amrex::Print() << " RHO*W RESIDUAL = " << rhow_residual << '\n';
      amrex::Print() << " RHO*E RESIDUAL = " << rhoE_residual << '\n';

      const int log_index = 0;
      if (log_index >= 0) {
        std::ostream& data_log2 = parent->DataLog(log_index);

        // Write the quantities at this time
        const int datwidth = 14;
        const int datprecision = 6;
        data_log2 << std::setw(datwidth) << time;
        data_log2 << std::setw(datwidth) << std::setprecision(datprecision)
                  << rho_mms_err;
        data_log2 << std::setw(datwidth) << std::setprecision(datprecision)
                  << u_mms_err;
        data_log2 << std::setw(datwidth) << std::setprecision(datprecision)
                  << v_mms_err;
        data_log2 << std::setw(datwidth) << std::setprecision(datprecision)
                  << w_mms_err;
        data_log2 << std::setw(datwidth) << std::setprecision(datprecision)
                  << p_mms_err;

        data_log2 << std::setw(datwidth) << std::setprecision(datprecision)
                  << rho_residual;
        data_log2 << std::setw(datwidth) << std::setprecision(datprecision)
                  << rhou_residual;
        data_log2 << std::setw(datwidth) << std::setprecision(datprecision)
                  << rhov_residual;
        data_log2 << std::setw(datwidth) << std::setprecision(datprecision)
                  << rhow_residual;
        data_log2 << std::setw(datwidth) << std::setprecision(datprecision)
                  << rhoE_residual;
        data_log2 << std::endl;
      }
    }
  }
}