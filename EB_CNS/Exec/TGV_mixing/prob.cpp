#include "prob.H"

extern "C" {
void amrex_probinit(const int* /*init*/, const int* /*name*/, const int* /*namelen*/,
                    const amrex::Real* problo, const amrex::Real* probhi)
{
  // Parse params
  {
    amrex::ParmParse pp("prob");
    pp.query("reynolds", CNS::h_prob_parm->reynolds);
    pp.query("mach", CNS::h_prob_parm->mach);
    pp.query("prandtl", CNS::h_prob_parm->prandtl);
    pp.query("convecting", CNS::h_prob_parm->convecting);
    pp.query("omega_x", CNS::h_prob_parm->omega_x);
    pp.query("omega_y", CNS::h_prob_parm->omega_y);
    pp.query("omega_z", CNS::h_prob_parm->omega_z);
  }

  // Define the length scale
  CNS::h_prob_parm->L = 1.0 / M_PI;
  // CNS::h_prob_parm->L_x = probhi[0] - problo[0];
  // CNS::h_prob_parm->L_y = probhi[1] - problo[1];
  // CNS::h_prob_parm->L_z = probhi[2] - problo[2];

  // Initial density, velocity, and material properties
  amrex::Real eint;
  amrex::Real cs;
  amrex::Real cp;
  amrex::Real massfrac[2] = {1.0, 0.0};
  auto eos = pele::physics::PhysicsType::eos();
  eos.PYT2RE(CNS::h_prob_parm->p0, massfrac, CNS::h_prob_parm->T0,
             CNS::h_prob_parm->rho0, eint);
  eos.RTY2Cs(CNS::h_prob_parm->rho0, CNS::h_prob_parm->T0, massfrac, cs);
  eos.TY2Cp(CNS::h_prob_parm->T0, massfrac, cp);
  CNS::h_prob_parm->v0 = CNS::h_prob_parm->mach * cs;

  auto& trans_parm = CNS::trans_parms.host_trans_parm();
  trans_parm.const_bulk_viscosity = 0.0;
  trans_parm.const_diffusivity = 0.0;
  trans_parm.const_viscosity = CNS::h_prob_parm->rho0 * CNS::h_prob_parm->v0 *
                               CNS::h_prob_parm->L / CNS::h_prob_parm->reynolds;
  trans_parm.const_conductivity =
    trans_parm.const_viscosity * cp / CNS::h_prob_parm->prandtl;
  CNS::trans_parms.sync_to_device();

  // Output IC
  amrex::Print() << "V0 = " << CNS::h_prob_parm->v0
                 << " p0 = " << CNS::h_prob_parm->p0
                 << " rho0 = " << CNS::h_prob_parm->rho0
                 << " T0 = " << CNS::h_prob_parm->T0 << std::endl;

  Gpu::copyAsync(Gpu::hostToDevice, CNS::h_prob_parm, CNS::h_prob_parm + 1,
                 CNS::d_prob_parm);
  Gpu::streamSynchronize();
}
}

void CNS::fill_ext_src(int i, int j, int k, amrex::Real /*time*/,
                       amrex::GeometryData const& /*geomdata*/,
                       amrex::Array4<const amrex::Real> const& state,
                       amrex::Array4<amrex::Real> const& ext_src,
                       Parm const& /*parm*/, ProbParm const& /*prob_parm*/)
{
}


void CNS::full_prob_post_timestep(int /*iteration*/)
{
  // Sum kinetic energy and incompressible enstrophy

  int finest_level = parent->finestLevel();
  amrex::Real time = state[State_Type].curTime();
  amrex::Real kinetic_e = 0.0;
  amrex::Real enstrophy = 0.0;
  amrex::Real k_sgs = 0.0;
  amrex::Real mixedness = 0.0;
  amrex::Real msgs = 0.0;

  if (level == 0) {
    if (amrex::ParallelDescriptor::IOProcessor()) {
      amrex::Print() << "... Calculating sum of quantities" << std::endl;
    }

    for (int lev = 0; lev <= finest_level; lev++) {
      CNS& cns_lev = getLevel(lev);

      // Sum kinetic energy
      amrex::MultiFab& S_new = cns_lev.get_new_data(State_Type);
      amrex::iMultiFab ifine_mask(cns_lev.grids, cns_lev.dmap, 1, 0);
      if (lev < parent->finestLevel()) {
        // mask out fine covered cells, do not sum
        ifine_mask = makeFineMask(cns_lev.grids, cns_lev.dmap, parent->boxArray(lev + 1), cns_lev.fine_ratio, 1, 0);
      } else {
        ifine_mask.setVal(1);
      }
      amrex::MultiFab volume(cns_lev.grids, cns_lev.dmap, 1, 0);
      cns_lev.geom.GetVolume(volume);
      auto mf = cns_lev.derive("magvort", time, 0);
      amrex::MultiFab magvort(cns_lev.grids, cns_lev.dmap, 1, 0);
      amrex::MultiFab::Copy(magvort, *mf, 0, 0, 1, 0);

      const auto geomdata = cns_lev.geom.data();
      auto const& sarrs = S_new.const_arrays();
      auto const& marrs = ifine_mask.const_arrays();
      auto const& volarr = volume.const_arrays();
      auto const& vortarrs = magvort.const_arrays();

      auto reduce_tuple = amrex::ParReduce(
        TypeList<ReduceOpSum, ReduceOpSum, ReduceOpSum, ReduceOpSum, ReduceOpSum>{}, 
        TypeList<Real, Real, Real, Real, Real>{}, S_new, IntVect(0),
        [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k) -> GpuTuple<Real, Real, Real, Real, Real> {
          const amrex::Real rho = sarrs[box_no](i, j, k, URHO);
          const amrex::Real rhoinv = 1.0 / rho;
          const amrex::Real ux = sarrs[box_no](i, j, k, UMX) * rhoinv;
          const amrex::Real uy = sarrs[box_no](i, j, k, UMY) * rhoinv;
          const amrex::Real uz = sarrs[box_no](i, j, k, UMZ) * rhoinv;
          const amrex::Real mask = amrex::Real(marrs[box_no](i, j, k));
          const amrex::Real vol = volarr[box_no](i, j, k);
          const amrex::Real mvt = vortarrs[box_no](i, j, k);
          const amrex::Real Y1 = sarrs[box_no](i, j, k, UFS) * rhoinv;

          amrex::Real ke = mask * vol * 0.5 * rho * (ux * ux + uy * uy + uz * uz);
          amrex::Real en = mask * vol * 0.5 * rho * (mvt * mvt); // 0.5 * rho * omega^2
          amrex::Real mxs = mask * vol * Y1 * (1.0 - Y1);

          amrex::Real ksgs = 0.0;
          amrex::Real mxsgs = 0.0;
#if NUM_FIELD > 1
          for (int nf = 1; nf <= NUM_FIELD; ++nf) {
            const amrex::Real rhof = sarrs[box_no](i, j, k, URHO + nf * NVAR);
            const amrex::Real rhofinv = 1.0 / rhof;
            const amrex::Real ufx = sarrs[box_no](i, j, k, UMX + nf * NVAR) * rhofinv;
            const amrex::Real ufy = sarrs[box_no](i, j, k, UMY + nf * NVAR) * rhofinv;
            const amrex::Real ufz = sarrs[box_no](i, j, k, UMZ + nf * NVAR) * rhofinv;
            const amrex::Real Yf1 = sarrs[box_no](i, j, k, UFS + nf * NVAR) * rhofinv;

            ksgs += mask * vol * 0.5 * rhof * ((ufx - ux) * (ufx - ux) 
                      + (ufy - uy) * (ufy - uy) + (ufz - uz) * (ufz - uz));
            mxsgs += mask * vol * Yf1 * (1.0 - Yf1);
          }
          ksgs /= amrex::Real(NUM_FIELD);
          mxsgs = mxsgs / amrex::Real(NUM_FIELD) - mxs;
#endif

          return {ke, en, ksgs, mxs, mxsgs};
        });
      kinetic_e += amrex::get<0>(reduce_tuple);
      enstrophy += amrex::get<1>(reduce_tuple);
      k_sgs += amrex::get<2>(reduce_tuple);
      mixedness += amrex::get<3>(reduce_tuple);
      msgs += amrex::get<4>(reduce_tuple);
    }

    // Reductions
    amrex::ParallelDescriptor::ReduceRealSum(
      &kinetic_e, 1, amrex::ParallelDescriptor::IOProcessorNumber());
    amrex::ParallelDescriptor::ReduceRealSum(
      &enstrophy, 1, amrex::ParallelDescriptor::IOProcessorNumber());
    amrex::ParallelDescriptor::ReduceRealSum(
      &k_sgs, 1, amrex::ParallelDescriptor::IOProcessorNumber());
    amrex::ParallelDescriptor::ReduceRealSum(
      &mixedness, 1, amrex::ParallelDescriptor::IOProcessorNumber());

    if (amrex::ParallelDescriptor::IOProcessor()) {
      amrex::Print() << "TIME = " << time << '\n'
                     << " SUM KE        = " << kinetic_e << " of which SGS = " << k_sgs << '\n'
                     << " SUM ENSTROPHY = " << enstrophy << '\n' // No SGS enstrophy in 1 point PDF
                     << " MIXEDNESS     = " << mixedness << " of which SGS = " << msgs << '\n';

      // Write the quantities at this time
      const int log_index = 0;
      std::ostream& data_log = parent->DataLog(log_index);
      const int datwidth = 14;
      const int datprecision = 6;
      data_log << std::setw(datwidth) << time;
      data_log << std::setw(datwidth) << std::setprecision(datprecision)
                << kinetic_e;
      data_log << std::setw(datwidth) << std::setprecision(datprecision)
                << enstrophy;
      data_log << std::setw(datwidth) << std::setprecision(datprecision)
                << k_sgs;
      data_log << std::setw(datwidth) << std::setprecision(datprecision)
                << mixedness;
      data_log << std::endl;      
    }
  }
}