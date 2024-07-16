#include "prob.H"

extern "C" {
void amrex_probinit(const int* /*init*/, const int* /*name*/, const int* /*namelen*/,
                    const amrex::Real* problo, const amrex::Real* probhi)
{
  // Parse params
  amrex::Real M = 0.1;
  {
    amrex::ParmParse pp("prob");
    pp.query("M", M);
    pp.query("T0", CNS::h_prob_parm->T0);
    pp.query("p0", CNS::h_prob_parm->p0);
  }

  // Premixed mass fraction
  CNS::h_prob_parm->Y[H2_ID] = 0.005833276514557218;
  CNS::h_prob_parm->Y[O2_ID] = 0.2315868836255521;
  CNS::h_prob_parm->Y[H2O_ID] = 3.0949523223470285e-05;
  CNS::h_prob_parm->Y[H_ID] = 9.41890468823083e-08;
  CNS::h_prob_parm->Y[O_ID] = 3.8043211459109995e-07;
  CNS::h_prob_parm->Y[OH_ID] = 1.590390948339689e-07;
  CNS::h_prob_parm->Y[HO2_ID] = 2.3400558812284596e-05;
  CNS::h_prob_parm->Y[H2O2_ID] = 2.074393015915917e-06;
  CNS::h_prob_parm->Y[N2_ID] = 0.762522781724583;

  // Define the length scale
  CNS::h_prob_parm->L = probhi[0] - problo[0];
  assert(CNS::h_prob_parm->L == probhi[1] - problo[1]);
  assert(CNS::h_prob_parm->L == probhi[2] - problo[2]);

  // Initialise velocity, density  
  amrex::Real cs, cp;
  auto eos = pele::physics::PhysicsType::eos();
  eos.PYT2R(CNS::h_prob_parm->p0, CNS::h_prob_parm->Y.begin(), CNS::h_prob_parm->T0,
            CNS::h_prob_parm->rho0);
  eos.RTY2Cs(CNS::h_prob_parm->rho0, CNS::h_prob_parm->T0,
             CNS::h_prob_parm->Y.begin(), cs);
  CNS::h_prob_parm->v0 = M * cs;
  eos.TY2Cp(CNS::h_prob_parm->T0, CNS::h_prob_parm->Y.begin(), cp);

  // Calculate Re and Pr
  amrex::Real mu, lam, dummy;
  bool get_mu = true, get_lam = true;
  auto trans = pele::physics::PhysicsType::transport();
  auto& trans_parm = CNS::trans_parms.host_trans_parm();
  trans.transport(false, get_mu, get_lam, false, false, CNS::h_prob_parm->T0,
                  CNS::h_prob_parm->rho0, CNS::h_prob_parm->Y.begin(), nullptr,
                  nullptr, mu, dummy, lam, &trans_parm);
  amrex::Real Re = CNS::h_prob_parm->rho0 * CNS::h_prob_parm->v0 *
                   CNS::h_prob_parm->L / mu;
  amrex::Real Pr = mu * cp / lam;
  amrex::Real t_ref = CNS::h_prob_parm->L / CNS::h_prob_parm->v0;

  // Print IC
  amrex::Print() << " ===== Initial conditions ===== \n"
                 << " v0 = " << CNS::h_prob_parm->v0 << '\n'
                 << " rho0 = " << CNS::h_prob_parm->rho0 << '\n'
                 << " Re = " << Re << '\n'
                 << " Pr = " << Pr << '\n'
                 << " t_ref = " << t_ref << '\n';

  Gpu::copy(Gpu::hostToDevice, CNS::h_prob_parm, CNS::h_prob_parm + 1,
            CNS::d_prob_parm);
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
  // Sum kinetic energy, incompressible enstrophy, mean temperature

  int finest_level = parent->finestLevel();
  amrex::Real time = state[State_Type].curTime();
  amrex::Real kinetic_e = 0.0;
  amrex::Real enstrophy = 0.0;
  amrex::Real mean_temp = 0.0;
  amrex::Real min_temp = 1.0e10;
  amrex::Real max_temp = 0.0;

  if (level == 0) {
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
      auto const& volarrs = volume.const_arrays();
      auto const& vortarrs = magvort.const_arrays();

      auto reduce_tuple = amrex::ParReduce(
        TypeList<ReduceOpSum, ReduceOpSum, ReduceOpSum, ReduceOpMin, ReduceOpMax>{}, 
        TypeList<Real, Real, Real, Real, Real>{}, S_new, IntVect(0),
        [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k) -> GpuTuple<Real, Real, Real, Real, Real> {
          const amrex::Real rho = sarrs[box_no](i, j, k, URHO);
          const amrex::Real rhoinv = 1.0 / rho;
          const amrex::Real mx = sarrs[box_no](i, j, k, UMX);
          const amrex::Real my = sarrs[box_no](i, j, k, UMY);
          const amrex::Real mz = sarrs[box_no](i, j, k, UMZ);
          const amrex::Real mask = amrex::Real(marrs[box_no](i, j, k));
          const amrex::Real vol = volarrs[box_no](i, j, k);
          const amrex::Real mvt = vortarrs[box_no](i, j, k);
          const amrex::Real ei = sarrs[box_no](i, j, k, UEDEN) * rhoinv -
                                 0.5 * (mx * mx + my * my + mz * mz) * rhoinv * rhoinv;
          amrex::Real Y[NUM_SPECIES];
          for (int n = 0; n < NUM_SPECIES; ++n)
            Y[n] = sarrs[box_no](i, j, k, UFS + n) * rhoinv;

          amrex::Real T;
          auto eos = pele::physics::PhysicsType::eos();
          eos.REY2T(rho, ei, Y, T);

          amrex::Real ke = mask * vol * 0.5 * rhoinv * (mx * mx + my * my + mz * mz);
          amrex::Real en = mask * vol * 0.5 * rho * (mvt * mvt); // 0.5 * rho * omega^2
          amrex::Real temp = mask * vol * T;
          return {ke, en, temp, T, T};
        });
      kinetic_e += amrex::get<0>(reduce_tuple);
      enstrophy += amrex::get<1>(reduce_tuple);
      mean_temp += amrex::get<2>(reduce_tuple);
      min_temp = amrex::min(min_temp, amrex::get<3>(reduce_tuple));
      max_temp = amrex::max(max_temp, amrex::get<4>(reduce_tuple));
    }

    // Reductions
    amrex::ParallelDescriptor::ReduceRealSum(
      &kinetic_e, 1, amrex::ParallelDescriptor::IOProcessorNumber());
    amrex::ParallelDescriptor::ReduceRealSum(
      &enstrophy, 1, amrex::ParallelDescriptor::IOProcessorNumber());
    amrex::ParallelDescriptor::ReduceRealSum(
      &mean_temp, 1, amrex::ParallelDescriptor::IOProcessorNumber());

    if (amrex::ParallelDescriptor::IOProcessor()) {
      amrex::Real vol = std::pow(CNS::h_prob_parm->L, 3);
      kinetic_e /= vol;
      enstrophy /= vol;
      mean_temp /= vol;

      amrex::Print() << "TIME = " << time << '\n'
                     << " SUM KE        = " << kinetic_e << '\n'
                     << " SUM ENSTROPHY = " << enstrophy << '\n'
                     << " TEMP          = " << min_temp << " .. " 
                                            << mean_temp << " .. " 
                                            << max_temp << '\n';

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
               << mean_temp;
      data_log << std::endl;
    }
  }
}