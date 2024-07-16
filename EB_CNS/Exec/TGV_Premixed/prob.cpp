#include "prob.H"
#include "file_io.H"

extern "C" {
void amrex_probinit(const int* /*init*/, const int* /*name*/, const int* /*namelen*/,
                    const amrex::Real* problo, const amrex::Real* probhi)
{
  // Parse params
  amrex::Real M = 0.1;
  std::string iname;
  bool binfmt = false;
  {
    amrex::ParmParse pp("prob");
    pp.query("M", M);
    pp.query("T0", CNS::h_prob_parm->T0);
    pp.query("p0", CNS::h_prob_parm->p0);
    pp.query("mode", CNS::h_prob_parm->mode);
    if (CNS::h_prob_parm->mode == 1) {
      pp.get("iname", iname);
      pp.get("inres", CNS::h_prob_parm->inres);
      pp.query("binfmt", binfmt);
    }
  }

  // Premixed mass fraction
  // CNS::h_prob_parm->Y[H2_ID] = 0.005833276514557218;
  // CNS::h_prob_parm->Y[O2_ID] = 0.2315868836255521;
  // CNS::h_prob_parm->Y[H2O_ID] = 3.0949523223470285e-05;
  // CNS::h_prob_parm->Y[H_ID] = 9.41890468823083e-08;
  // CNS::h_prob_parm->Y[O_ID] = 3.8043211459109995e-07;
  // CNS::h_prob_parm->Y[OH_ID] = 1.590390948339689e-07;
  // CNS::h_prob_parm->Y[HO2_ID] = 2.3400558812284596e-05;
  // CNS::h_prob_parm->Y[H2O2_ID] = 2.074393015915917e-06;
  // CNS::h_prob_parm->Y[N2_ID] = 0.762522781724583;
  // CNS::h_prob_parm->Y[H2_ID] = 0.0058188147027386705;
  // CNS::h_prob_parm->Y[O2_ID] = 0.23146147387598764;
  // CNS::h_prob_parm->Y[H2O_ID] = 0.00014896541997113147;
  // CNS::h_prob_parm->Y[H_ID] = 6.646291735036124e-07;
  // CNS::h_prob_parm->Y[O_ID] = 2.333958270384001e-06;
  // CNS::h_prob_parm->Y[OH_ID] = 1.4468188194152512e-06;
  // CNS::h_prob_parm->Y[HO2_ID] = 3.94531289793572e-05;
  // CNS::h_prob_parm->Y[H2O2_ID] = 4.065741477139064e-06;
  // CNS::h_prob_parm->Y[N2_ID] = 0.7625227817245828;
  // CNS::h_prob_parm->Y[H2_ID] = 0.02849334158827717;
  // CNS::h_prob_parm->Y[O2_ID] = 0.22607628141952604;
  // CNS::h_prob_parm->Y[H2O_ID] = 0.0002165116676253062;
  // CNS::h_prob_parm->Y[H_ID] = 1.0509296964174552e-06;
  // CNS::h_prob_parm->Y[O_ID] = 9.443135255685144e-07;
  // CNS::h_prob_parm->Y[OH_ID] = 5.710268033998482e-07;
  // CNS::h_prob_parm->Y[HO2_ID] = 5.1010518050022254e-05;
  // CNS::h_prob_parm->Y[H2O2_ID] = 3.66830350703223e-05;
  // CNS::h_prob_parm->Y[N2_ID] = 0.7451236055014252;
  // FakeReaction
  // CNS::h_prob_parm->Y[H2_ID] = 0.02;
  // CNS::h_prob_parm->Y[O2_ID] = 0.98;
  // CNS::h_prob_parm->Y[H2O_ID] = 0.0;
  // FakeReaction2
  CNS::h_prob_parm->Y[F1_ID] = 0.5;
  CNS::h_prob_parm->Y[F2_ID] = 0.5;

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

  // HIT
  bool restart = false;
  {
    amrex::ParmParse pp("amr");
    std::string start_file;
    pp.query("restart", start_file);
    if (!start_file.empty()) 
      restart = true;
  }

  if (!restart && CNS::h_prob_parm->mode == 1) {
#ifdef AMREX_USE_FLOAT
    amrex::Abort("HIT cannot run in single precision at the moment.");
#else
// #if NUM_FIELD == 0
//     int num_entry = 6; // x, y, z, u, v, w
// #else
    int num_entry = 7; // ..., k_sgs
// #endif
    const size_t nx = CNS::h_prob_parm->inres;
    const size_t ny = CNS::h_prob_parm->inres;
    const size_t nz = CNS::h_prob_parm->inres;
    amrex::Vector<amrex::Real> data(nx * ny * nz * num_entry);
    if (binfmt) {
      read_binary(iname, nx, ny, nz, num_entry, data);
    } else {
      read_csv(iname, nx, ny, nz, data);
    }

    // Extract position and velocities
    amrex::Vector<amrex::Real> h_xinput(nx * ny * nz);
    CNS::prob_parm_host->h_uinput.resize(nx * ny * nz);
    CNS::prob_parm_host->h_vinput.resize(nx * ny * nz);
    CNS::prob_parm_host->h_winput.resize(nx * ny * nz);
// #if NUM_FIELD > 0
    CNS::prob_parm_host->h_kinput.resize(nx * ny * nz);
// #endif
    for (long i = 0; i < h_xinput.size(); i++) {
      h_xinput[i] = data[0 + i * num_entry];
      CNS::prob_parm_host->h_uinput[i] = data[3 + i * num_entry] * CNS::h_prob_parm->v0;
      CNS::prob_parm_host->h_vinput[i] = data[4 + i * num_entry] * CNS::h_prob_parm->v0;
      CNS::prob_parm_host->h_winput[i] = data[5 + i * num_entry] * CNS::h_prob_parm->v0;
// #if NUM_FIELD > 0
      CNS::prob_parm_host->h_kinput[i] = data[6 + i * num_entry] * CNS::h_prob_parm->v0 * CNS::h_prob_parm->v0;
      CNS::h_prob_parm->ksgs_avg += CNS::prob_parm_host->h_kinput[i] / h_xinput.size();
// #endif
    }

    // Get the xarray table and the differences.
    amrex::Vector<amrex::Real> h_xarray(nx);
    for (long i = 0; i < h_xarray.size(); i++) {
      h_xarray[i] = h_xinput[i];
    }
    amrex::Vector<amrex::Real> h_xdiff(nx);
    std::adjacent_difference(h_xarray.begin(), h_xarray.end(), h_xdiff.begin());
    h_xdiff[0] = h_xdiff[1];

    // Make sure the search array is increasing
    if (!std::is_sorted(h_xarray.begin(), h_xarray.end())) {
      amrex::Abort("Error: non ascending x-coordinate array.");
    }

    // Get pointer to the data
    CNS::prob_parm_host->uinput.resize(CNS::prob_parm_host->h_uinput.size());
    CNS::prob_parm_host->vinput.resize(CNS::prob_parm_host->h_vinput.size());
    CNS::prob_parm_host->winput.resize(CNS::prob_parm_host->h_winput.size());
// #if NUM_FIELD > 0
    CNS::prob_parm_host->kinput.resize(CNS::prob_parm_host->h_kinput.size());
// #endif
    amrex::Gpu::copyAsync(
      amrex::Gpu::hostToDevice, CNS::prob_parm_host->h_uinput.begin(),
      CNS::prob_parm_host->h_uinput.end(), CNS::prob_parm_host->uinput.begin());
    amrex::Gpu::copyAsync(
      amrex::Gpu::hostToDevice, CNS::prob_parm_host->h_vinput.begin(),
      CNS::prob_parm_host->h_vinput.end(), CNS::prob_parm_host->vinput.begin());
    amrex::Gpu::copyAsync(
      amrex::Gpu::hostToDevice, CNS::prob_parm_host->h_winput.begin(),
      CNS::prob_parm_host->h_winput.end(), CNS::prob_parm_host->winput.begin());
// #if NUM_FIELD > 0
    amrex::Gpu::copyAsync(
      amrex::Gpu::hostToDevice, CNS::prob_parm_host->h_kinput.begin(),
      CNS::prob_parm_host->h_kinput.end(), CNS::prob_parm_host->kinput.begin());
// #endif    
    CNS::h_prob_parm->d_uinput = CNS::prob_parm_host->uinput.data();
    CNS::h_prob_parm->d_vinput = CNS::prob_parm_host->vinput.data();
    CNS::h_prob_parm->d_winput = CNS::prob_parm_host->winput.data();
// #if NUM_FIELD > 0
    CNS::h_prob_parm->d_kinput = CNS::prob_parm_host->kinput.data();
// #endif
#endif
  }

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
        TypeList<ReduceOpSum, ReduceOpSum, ReduceOpSum, ReduceOpMax, ReduceOpMin>{},
        TypeList<Real, Real, Real, Real, Real>{}, S_new, IntVect(0),
        [=] AMREX_GPU_DEVICE(int box_no, int i, int j,
                             int k) -> GpuTuple<Real, Real, Real, Real, Real> {
          if (marrs[box_no](i, j, k) == 0) {
            return {0.0, 0.0, 0.0, 1.0e10, 0.0};
          }

          const amrex::Real vol = volarrs[box_no](i, j, k);
          const amrex::Real mvt = vortarrs[box_no](i, j, k);

          amrex::Real rho = sarrs[box_no](i, j, k, URHO);
          amrex::Real rhoinv = 1.0 / rho;
          amrex::Real mx = sarrs[box_no](i, j, k, UMX);
          amrex::Real my = sarrs[box_no](i, j, k, UMY);
          amrex::Real mz = sarrs[box_no](i, j, k, UMZ);
          amrex::Real ei = sarrs[box_no](i, j, k, UEDEN) * rhoinv -
                           0.5 * (mx * mx + my * my + mz * mz) * rhoinv * rhoinv;
          amrex::Real Y[NUM_SPECIES];
          for (int n = 0; n < NUM_SPECIES; ++n)
            Y[n] = sarrs[box_no](i, j, k, UFS + n) * rhoinv;

          amrex::Real T;
          auto eos = pele::physics::PhysicsType::eos();
          eos.REY2T(rho, ei, Y, T);

          amrex::Real ke = vol * 0.5 * rhoinv * (mx * mx + my * my + mz * mz);
          amrex::Real en = vol * 0.5 * rho * (mvt * mvt); // 0.5 * rho * omega^2
          amrex::Real temp = vol * T;

          amrex::Real maxT = T;
          amrex::Real minT = T;
#if NUM_FIELD > 0
          amrex::Real ksgs = 0.0;
          for (int nf = 1; nf <= NUM_FIELD; ++nf) {
            // Temperature of each field
            amrex::Real rhof = sarrs[box_no](i, j, k, nf * NVAR + URHO);
            amrex::Real rhofinv = 1.0 / rhof;
            amrex::Real mxf = sarrs[box_no](i, j, k, nf * NVAR + UMX);
            amrex::Real myf = sarrs[box_no](i, j, k, nf * NVAR + UMY);
            amrex::Real mzf = sarrs[box_no](i, j, k, nf * NVAR + UMZ);
            ei = sarrs[box_no](i, j, k, nf * NVAR + UEDEN) * rhofinv -
                 0.5 * (mxf * mxf + myf * myf + mzf * mzf) * rhofinv * rhofinv;
            for (int n = 0; n < NUM_SPECIES; ++n)
              Y[n] = sarrs[box_no](i, j, k, nf * NVAR + UFS + n) * rhofinv;
            eos.REY2T(rhof, ei, Y, T);
            maxT = amrex::max(maxT, T);
            minT = amrex::min(minT, T);

            ksgs += vol * 0.5 * rho * (pow(mxf * rhofinv - mx * rhoinv, 2) 
                                     + pow(myf * rhofinv - my * rhoinv, 2) 
                                     + pow(mzf * rhofinv - mz * rhoinv, 2));
          }
          ksgs /= amrex::Real(NUM_FIELD);
          ke += ksgs;
#endif
          return {ke, en, temp, maxT, minT};
        });
      kinetic_e += amrex::get<0>(reduce_tuple);
      enstrophy += amrex::get<1>(reduce_tuple);
      mean_temp += amrex::get<2>(reduce_tuple);
      max_temp = amrex::max(max_temp, amrex::get<3>(reduce_tuple));
      min_temp = amrex::min(min_temp, amrex::get<4>(reduce_tuple));
    }

    // Reductions
    amrex::ParallelDescriptor::ReduceRealSum(
      &kinetic_e, 1, amrex::ParallelDescriptor::IOProcessorNumber());
    amrex::ParallelDescriptor::ReduceRealSum(
      &enstrophy, 1, amrex::ParallelDescriptor::IOProcessorNumber());
    amrex::ParallelDescriptor::ReduceRealSum(
      &mean_temp, 1, amrex::ParallelDescriptor::IOProcessorNumber());
    amrex::ParallelDescriptor::ReduceRealMax(
      &max_temp, 1, amrex::ParallelDescriptor::IOProcessorNumber());
    amrex::ParallelDescriptor::ReduceRealMin(
      &min_temp, 1, amrex::ParallelDescriptor::IOProcessorNumber());

    if (amrex::ParallelDescriptor::IOProcessor()) {
      amrex::Real vol = std::pow(CNS::h_prob_parm->L, 3);
      kinetic_e /= vol;
      enstrophy /= vol;
      mean_temp /= vol;

      amrex::Print() << "TIME = " << time << '\n'
                     << " SUM KE        = " << kinetic_e << '\n'
                     << " SUM ENSTROPHY = " << enstrophy << '\n'
                     << " TEMP          = " << min_temp << ".." 
                                            << mean_temp << ".." 
                                            << max_temp << '\n';

      // Write the quantities at this time
      const int log_index = 0;
      std::ostream& data_log = parent->DataLog(log_index);
      const int datwidth = 14;
      const int datprecision = 6;
      data_log << std::setw(datwidth) << time;
      data_log << std::setw(datwidth) << std::setprecision(datprecision) << kinetic_e;
      data_log << std::setw(datwidth) << std::setprecision(datprecision) << enstrophy;
      data_log << std::setw(datwidth) << std::setprecision(datprecision) << mean_temp;
      data_log << std::setw(datwidth) << std::setprecision(datprecision) << min_temp;
      data_log << std::setw(datwidth) << std::setprecision(datprecision) << max_temp;
      data_log << std::endl;
    }
  }
}