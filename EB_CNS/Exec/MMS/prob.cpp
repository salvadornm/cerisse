#include "prob.H"

extern "C" {
void amrex_probinit(const int* /*init*/, const int* /*name*/, const int* /*namelen*/,
                    const amrex::Real* problo, const amrex::Real* probhi)
{
  // Parse params
  std::string masa_solution_name = "navierstokes_3d_compressible";
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
    pp.query("masa_solution_name", masa_solution_name);
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
  trans_parm.const_bulk_viscosity = 0.0;
  trans_parm.const_diffusivity = 0.0;
  trans_parm.const_viscosity = CNS::h_prob_parm->rho0 * CNS::h_prob_parm->u0 *
                               CNS::h_prob_parm->L / CNS::h_prob_parm->reynolds;
  trans_parm.const_conductivity =
    trans_parm.const_viscosity * cp / CNS::h_prob_parm->prandtl;
  CNS::trans_parms.sync_to_device();

  Gpu::copyAsync(Gpu::hostToDevice, CNS::h_prob_parm, CNS::h_prob_parm + 1,
                 CNS::d_prob_parm);
  Gpu::streamSynchronize();

  // MASA parameters for the following functions
  // rho = rho_0 + rho_x * cos(a_rhox * PI * x / L) * rho_y * cos(a_rhoy * PI * y / L) * rho_z * cos(a_rhoz * PI * z / L)
  // u = u_0 + u_x * cos(a_ux * PI * x / L) * u_y * cos(a_uy * PI * y / L) * u_z * cos(a_uz * PI * z / L)
  // v = v_0 + v_x * cos(a_vx * PI * x / L) * v_y * cos(a_vy * PI * y / L) * v_z * cos(a_vz * PI * z / L)
  // w = w_0 + w_x * cos(a_wx * PI * x / L) * w_y * cos(a_wy * PI * y / L) * w_z * cos(a_wz * PI * z / L)
  // p = p_0 + p_x * cos(a_px * PI * x / L) * p_y * cos(a_py * PI * y / L) * p_z * cos(a_pz * PI * z / L)
  masa_init("mms", masa_solution_name.c_str());
  masa_set_param("L", CNS::h_prob_parm->L);
  masa_set_param("R",
                 pele::physics::Constants::RU / pele::physics::Constants::AIRMW);
  masa_set_param("k", trans_parm.const_conductivity);
  masa_set_param("Gamma", eos.gamma);
  masa_set_param("mu", trans_parm.const_viscosity);
  // masa_set_param("mu_bulk", trans_parm.const_bulk_viscosity);
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

  if (amrex::ParallelDescriptor::IOProcessor()) {
    masa_display_param();
  }
  masa_sanity_check();
}
}

void CNS::fill_ext_src(int i, int j, int k, amrex::Real time,
                       amrex::GeometryData const& geomdata,
                       amrex::Array4<const amrex::Real> const& state,
                       amrex::Array4<amrex::Real> const& ext_src,
                       Parm const& /*parm*/, ProbParm const& /*prob_parm*/)
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

void CNS::prob_post_coarsetimestep(amrex::Real time)
{
  if (level > 0) return;

  if (amrex::ParallelDescriptor::IOProcessor())
    amrex::Print() << "Calculating MMS errors" << std::endl;

  // Calculate errors
  amrex::MultiFab& S = get_new_data(State_Type);
  amrex::iMultiFab ifine_mask(grids, dmap, 1, 0);
  if (level < parent->finestLevel()) {
    ifine_mask = makeFineMask(grids, dmap, parent->boxArray(level + 1), fine_ratio, 1, 0);
  } else {
    ifine_mask.setVal(1);
  }    
  const auto geomdata = geom.data();
  auto const& sarrs = S.const_arrays();
  auto const& marrs = ifine_mask.const_arrays();

  GpuTuple<Real, Real, Real, Real, Real> errors = amrex::ParReduce(
    TypeList<ReduceOpSum, ReduceOpSum, ReduceOpSum, ReduceOpSum, ReduceOpSum>{},
    TypeList<Real, Real, Real, Real, Real>{}, S, IntVect(0),
    [=] AMREX_GPU_DEVICE(int box_no, int i, int j,
                         int k) -> GpuTuple<Real, Real, Real, Real, Real> {
      const amrex::Real* prob_lo = geomdata.ProbLo();
      const amrex::Real* dx = geomdata.CellSize();
      const amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];
      const amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1];
      const amrex::Real z = prob_lo[2] + (k + 0.5) * dx[2];

      // Get exact solution
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
      amrex::Real Y[NUM_SPECIES] = {0.0};
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

      return {rho_err, u_err, v_err, w_err, p_err};
    });

  amrex::Real rho_mms_err = amrex::get<0>(errors);
  amrex::Real u_mms_err = amrex::get<1>(errors);
  amrex::Real v_mms_err = amrex::get<2>(errors);
  amrex::Real w_mms_err = amrex::get<3>(errors);
  amrex::Real p_mms_err = amrex::get<4>(errors);

  if (amrex::ParallelDescriptor::IOProcessor()) {
    amrex::Real time = state[State_Type].curTime();
    amrex::Print() << "TIME= " << time << " RHO MMS ERROR = " << rho_mms_err << '\n';
    amrex::Print() << "TIME= " << time << " U MMS ERROR   = " << u_mms_err << '\n';
    amrex::Print() << "TIME= " << time << " V MMS ERROR   = " << v_mms_err << '\n';
    amrex::Print() << "TIME= " << time << " W MMS ERROR   = " << w_mms_err << '\n';
    amrex::Print() << "TIME= " << time << " P MMS ERROR   = " << p_mms_err << '\n';

    // Write the errors to the datalog
    // std::string data_log_name = "mms.log";
    // Find data_log index
    // int log_index = -1;
    // for (int ii = 0; ii < parent->NumDataLogs(); ii++) {
    //   if (parent->DataLogName(ii) == data_log_name) {
    //     log_index = ii;
    //   }
    // }
    int log_index = 0;
    // Write data to the log file
    if (log_index >= 0) {
      std::ostream& data_log2 = parent->DataLog(log_index);
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
      data_log2 << std::endl;
    }
  }
}