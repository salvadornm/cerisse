#include "CNS.H"
#include "central_scheme.H"
#include "diffusion.H"
#include "hydro.H"
#include "recon.H"
#include "wall_model.H"

using namespace amrex;

void CNS::compute_dSdt_box(Box const& bx, Array4<const Real>& sarr,
                           Array4<Real>& dsdt,
                           std::array<FArrayBox*, AMREX_SPACEDIM> const& flxfab,
                           Array4<const Real>& shock_sensor)
{
  BL_PROFILE("CNS::compute_dSdt_box()");

#if NUM_FIELD > 0
  // auto const& geomdata = geom.data();
  // const Real x = geomdata.ProbLo(0) + (bx.bigEnd(0) + 0.5) * geomdata.CellSize(0);
  // const bool update_fields = (x > -1.0) && (x < 5.0); // Do SF or just the mean? TODO: make this a tagging
  const bool update_fields = true;
#else
  const bool update_fields = false;
#endif
  const int ncomp = update_fields ? UFA : NVAR; // we don't want to change time averages
  const Box& bxg2 = amrex::grow(bx, 2);
  const Box& bxg3 = amrex::grow(bx, 3);
  const auto dx = geom.CellSizeArray();
  const auto dxinv = geom.InvCellSizeArray();
  const bool do_diffusion = do_visc || do_les || buffer_box.ok();

  // Prepare FABs to store data
  FArrayBox qfab(bxg3, NPRIM, The_Async_Arena()); // Primitive variables
  FArrayBox wfab(bxg3, NCHAR, The_Async_Arena()); // FABs for reconstruction
  FArrayBox wlfab(bxg2, NCHAR, The_Async_Arena());
  FArrayBox wrfab(bxg2, NCHAR, The_Async_Arena());
  FArrayBox coefsfab; // Diffusion coefficients
  if (do_diffusion) { coefsfab.resize(bxg2, NCOEF, The_Async_Arena()); }
  auto const& q = qfab.array();
  auto const& w = wfab.array();
  auto const& wl = wlfab.array();
  auto const& wr = wrfab.array();
  auto const& coefs = coefsfab.array();

  // Store viscous fluxes separately in V/VSPDF
  std::array<FArrayBox, amrex::SpaceDim> vfluxfab;
  const bool store_in_vflux = do_diffusion && do_vpdf && update_fields;
  if (store_in_vflux) {
    for (int dir = 0; dir < amrex::SpaceDim; ++dir) {
      vfluxfab[dir].resize(amrex::surroundingNodes(bx, dir), NVAR,
                           The_Async_Arena());
      vfluxfab[dir].setVal<RunOn::Device>(0.0);
    }
  }

  // Advance
#if NUM_FIELD > 0
  const int nf_start = update_fields ? 1 : 0;
  const int nf_end = update_fields ? NUM_FIELD : 0;
  for (int nf = nf_start; nf <= nf_end; ++nf)
#else
  const int nf = 0;
#endif
  {
    // Prepare primitive variables
    amrex::ParallelFor(bxg3, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
      cns_ctoprim(i, j, k, nf * NVAR, sarr, q);
    });

    // Prepare transport coefs
    if (do_diffusion) {
      auto const& qar_yin = qfab.const_array(QFS);
      auto const& qar_Tin = qfab.const_array(QTEMP);
      auto const& qar_rhoin = qfab.const_array(QRHO);
      auto const& mu = coefsfab.array(CMU);      // dynamic viscosity
      auto const& xi = coefsfab.array(CXI);      // bulk viscosity
      auto const& lambda = coefsfab.array(CLAM); // thermal conductivity
      auto const& rhoD = coefsfab.array(CRHOD);  // species diffusivity (times rho)

      // Physical transport coefs
      if (do_visc) {
        BL_PROFILE("PelePhysics::get_transport_coeffs()");
        Array4<Real> chi; // dummy Soret effect coef
        auto const* ltransparm = trans_parms.device_trans_parm();
        amrex::launch(bxg2, [=] AMREX_GPU_DEVICE(Box const& tbx) {
          auto trans = pele::physics::PhysicsType::transport();
          trans.get_transport_coeffs(tbx, qar_yin, qar_Tin, qar_rhoin, rhoD, chi, mu,
                                     xi, lambda, ltransparm);
        });
      }

      // LES diffucsion coefs
      if (do_les) {
        BL_PROFILE("CNS::LES_transport_coeffs");
        Real delta = std::pow(AMREX_D_TERM(dx[0], *dx[1], *dx[2]),
                              Real(1.0) / Real(amrex::SpaceDim)); // LES filter width
        amrex::ParallelFor(bxg2, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
          Real mu_T, xi_T;
          les_model->mu_T_cc(i, j, k, qar_rhoin, dxinv, delta, Cs, mu_T);
          // les_model->xi_T_cc(i, j, k, qar_rhoin, dxinv, delta, C_I, xi_T);
          xi_T = 0.0;
          // TODO: this divu is calculated at cell centre, not face centre
          Real Y[NUM_SPECIES], cp;
          for (int n = 0; n < NUM_SPECIES; ++n) { Y[n] = qar_yin(i, j, k, n); }
          auto eos = pele::physics::PhysicsType::eos();
          eos.RTY2Cp(qar_rhoin(i, j, k), qar_Tin(i, j, k), Y, cp);

          mu(i, j, k) += mu_T;
          xi(i, j, k) += xi_T;
          lambda(i, j, k) += cp * mu_T / Pr_T;
          for (int ns = 0; ns < NUM_SPECIES; ++ns) {
            rhoD(i, j, k, ns) += mu_T / Sc_T;
          }
        });
      }

      // Buffer region
      if (buffer_box.ok()) {
        const auto problo = geom.ProbLo();
        amrex::ParallelFor(bxg2, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
          RealVect pos{AMREX_D_DECL((i + 0.5) * dx[0] + problo[0],
                                    (j + 0.5) * dx[1] + problo[1],
                                    (k + 0.5) * dx[2] + problo[2])};
          if (buffer_box.contains(pos)) {
            mu(i, j, k) += 0.1;
            xi(i, j, k) += 0.1;
            lambda(i, j, k) += 1e5;
          }
        });
      }
    } // if do_diffusion

    int num_scheme_switch = 0;
    if (use_hybrid_scheme) {
      ReduceOps<ReduceOpMax> reduce_op;
      ReduceData<Real> reduce_data(reduce_op);
      using ReduceTuple = typename decltype(reduce_data)::Type;
      reduce_op.eval(bxg2, reduce_data,
                     [=] AMREX_GPU_DEVICE(int i, int j, int k) -> ReduceTuple {
                       return {shock_sensor(i, j, k)};
                     });
      Real max_shock_sensor = amrex::get<0>(reduce_data.value());

      if (max_shock_sensor < 0.7) {
        num_scheme_switch = 1; // all smooth, use central scheme
                               // else, use shock-capturing scheme
      }
    }

    // Compute fluxes for each space direction
    for (int dir = 0; dir < amrex::SpaceDim; ++dir) {
      const Box& flxbx = amrex::surroundingNodes(bx, dir);
      auto const& flx = flxfab[dir]->array(nf * NVAR);

      // Hydro/hyperbolic fluxes
      if (do_hydro) {
        if (num_scheme_switch == 0) {
          // Shock-capturing scheme
          // 1. Convert primitive to characteristic at cell centre
          const Box& charbox = amrex::grow(bx, dir, 3);
          amrex::ParallelFor(charbox, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            cns_ctochar(i, j, k, dir, q, w, char_sys);
          });
          // 2. FD interpolation to cell face
          const Box& reconbox = amrex::grow(bx, dir, 1);
          amrex::ParallelFor(TypeList<CompileTimeOptions<1, 2, 3, 4, 5, 6>>{},
                             {recon_scheme}, reconbox, NCHAR,
                             [=] AMREX_GPU_DEVICE(int i, int j, int k, int n,
                                                  auto captured_recon_scheme) {
                               cns_recon<captured_recon_scheme>(i, j, k, n, dir, w,
                                                                wl, wr, plm_theta);
                             });
          // 3. Solve Riemann problem for fluxes at cell face
          amrex::ParallelFor(flxbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            cns_riemann(i, j, k, dir, flx, q, wl, wr, char_sys, recon_char_var);

            // bool do_high_order_diff =
            //   (shock_sensor(i, j, k) < 0.95) &&
            //   (shock_sensor(IntVect(AMREX_D_DECL(i, j, k)) -
            //                 IntVect::TheDimensionVector(dir)) < 0.95);
            // if (do_high_order_diff) { cns_afd_correction(i, j, k, dir, q, flx); }
          });
        } else {
          // Central scheme
          amrex::ParallelFor(flxbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            cns_KEEP4(i, j, k, dir, flx, q);
          });
        }
      }
     
      // Viscous fluxes
      if (do_diffusion) {
        auto const& vflx = store_in_vflux ? vfluxfab[dir].array() : flx;
        
        amrex::ParallelFor(flxbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
          const IntVect iv(AMREX_D_DECL(i, j, k));

          Real flx_tmp[NVAR] = {0.0};
          cns_diff(iv, dir, q, coefs, dxinv, flx_tmp);
          
          // Wall model for regular solid boundaries (modifies flx_tmp)
          const bool lo_is_wall = phys_bc.lo(dir) == 5;
          const bool hi_is_wall = phys_bc.hi(dir) == 5;
          const int domlo = geom.Domain().smallEnd(dir);
          const int domhi = geom.Domain().bigEnd(dir);
          const auto problo = geom.ProbLo();
          // const Real x = (i + 0.5) * dx[0] + problo[0];
          if (eb_wall_model /*&& x < -1.0*/) { // in line with advance_box_eb            
            if ((iv[dir] == domhi + 1 && hi_is_wall) ||
                (iv[dir] == domlo && lo_is_wall)) {
              const auto iv1 = iv[dir] == domlo ? iv : iv - IntVect::TheDimensionVector(dir);
              
              // setup some image point values
              // get density and transport coefficients
              Real rho, mu, lam, p = q(iv1, QPRES),
                                 Y[NUM_SPECIES]; // assume dp/dn, dY/dn ~ 0
              for (int n = 0; n < NUM_SPECIES; ++n) { Y[n] = q(iv1, QFS + n); }
              // if (eb_isothermal) {
              //   auto trans = pele::physics::PhysicsType::transport();
              //   auto const* ltransparm = CNS::trans_parms.device_trans_parm();
              //   Real xi_unused;

              //   auto eos = pele::physics::PhysicsType::eos();
              //   eos.PYT2R(p, Y, eb_wall_temp, rho); // assume dp/dn ~ 0

              //   trans.transport(true, true, true, false, false, eb_wall_temp, rho, Y,
              //                   nullptr, nullptr, mu, xi_unused, lam, ltransparm);
              // } else 
              {
                // Adiabatic wall dT/dn ~ 0
                rho = q(iv1, QRHO);
                mu = coefs(iv1, CMU);
                lam = coefs(iv1, CLAM);
              }

              // Check y+
              const Real local_dx = dx[0];  // dx == dy == dz
              const Real tau =
                std::sqrt(flx_tmp[UMX] * flx_tmp[UMX] + flx_tmp[UMY] * flx_tmp[UMY] +
                          flx_tmp[UMZ] * flx_tmp[UMZ]);
              const Real utau = std::sqrt(tau / rho);
              const Real yplus = rho * utau / mu * local_dx;              

              // Tangent vectors (t1.n = 0, t2.t1 = 0, u.t2 = 0)
              const auto iv2 = iv[dir] == domlo ? iv + IntVect::TheDimensionVector(dir)
                                                : iv - 2 * IntVect::TheDimensionVector(dir);
              // const auto iv2 = iv + lohi * IntVect::TheDimensionVector(dir);
              Real u2[3] = {q(iv2, QU), q(iv2, QV), q(iv2, QW)};
              u2[dir] = 0.0;
              const Real u_parallel = std::sqrt(u2[0] * u2[0] + u2[1] * u2[1] + u2[2] * u2[2]);
              const Real t1[3] = {u2[0] / u_parallel, u2[1] / u_parallel, u2[2] / u_parallel};
              const Real ts = (u2[0] > 0.0) ? 1.0 : -1.0;
              const Real T2 = q(iv2, QTEMP);
                            
              // call wall_model.parallel_wall_stress
              if (yplus > 11.0 && u_parallel > 100.0 && !(T2 < 90.0) &&
                  !(T2 > 4000.0) && rho > 0.0) {
                // LawOfTheWall wm;
                EquilibriumODE wm;                
                Real T_wall = eb_isothermal ? eb_wall_temp : -1.0;
                Real h = 1.5 * local_dx;
                Real tauw, qw;
                wm.parallel_wall_stress(u_parallel, T2, rho, Y, h, mu, lam, T_wall,
                                        tauw, qw);
                if (!isnan(tauw) && !isnan(qw)) {
                  const Real lohi = iv[dir] == domlo ? -1 : 1;
                  flx_tmp[UMX] = ts * t1[0] * tauw * lohi;
                  flx_tmp[UMY] = ts * t1[1] * tauw * lohi;
                  flx_tmp[UMZ] = ts * t1[2] * tauw * lohi;
                  flx_tmp[UEDEN] = qw * lohi; // TODO: no wall model heat flux for now
                }
              }
            }
          }

          // Add flx_tmp to vflx
          for (int n = 0; n < NVAR; ++n) { vflx(iv, n) += flx_tmp[n]; }
        });
      }
    } // for dir
  } // for fields

#if NUM_FIELD > 0
  // Average viscous fluxes for V/VSPDF and add to flx
  if (store_in_vflux) {
    for (int dir = 0; dir < amrex::SpaceDim; ++dir) {
      const Box& flxbx = amrex::surroundingNodes(bx, dir);
      auto const& flx = flxfab[dir]->array();
      auto const& vflx = vfluxfab[dir].array();

      amrex::ParallelFor(
        flxbx, NUM_FIELD, [=] AMREX_GPU_DEVICE(int i, int j, int k, int nfm1) {
          const Real invNF = Real(1.0) / Real(NUM_FIELD);
          const int nf = 1 + nfm1;
          AMREX_D_TERM(flx(i, j, k, nf * NVAR + UMX) += invNF * vflx(i, j, k, UMX);
                       , flx(i, j, k, nf * NVAR + UMY) += invNF * vflx(i, j, k, UMY);
                       , flx(i, j, k, nf * NVAR + UMZ) += invNF * vflx(i, j, k, UMZ);)
          flx(i, j, k, nf * NVAR + UEDEN) += invNF * vflx(i, j, k, UEDEN);
          for (int n = 0; n < NUM_SPECIES; ++n) {
            flx(i, j, k, nf * NVAR + UFS + n) += invNF * vflx(i, j, k, UFS + n);
          }
        });
    } // for dir
  }
#endif

  // Compute flux divergence
  amrex::ParallelFor(bx, ncomp,
                     [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
                       cns_div(i, j, k, n, dsdt,
                               AMREX_D_DECL(flxfab[0]->array(), flxfab[1]->array(),
                                            flxfab[2]->array()),
                               dxinv);
                     });

#if NUM_FIELD > 0
{
  // BL_PROFILE("New_field_averaging");
  // Copy mean field to all fields if not doing SF, average to mean if doing SF
  for (int dir = 0; dir < amrex::SpaceDim; ++dir) {
    const Box& flxbx = amrex::surroundingNodes(bx, dir);
    auto const& flx = flxfab[dir]->array();

    amrex::ParallelFor(
      flxbx, NUM_FIELD, [=] AMREX_GPU_DEVICE(int i, int j, int k, int nfm1) {
        const int nf = 1 + nfm1;
        const Real fac = update_fields ? Real(1.0) / Real(NUM_FIELD) : 1.0;
        const int from_id = update_fields ? nf * NVAR : 0;
        const int to_id = update_fields ? 0 : nf * NVAR;

        flx(i, j, k, to_id + URHO) += fac * flx(i, j, k, from_id + URHO);
        AMREX_D_TERM(
          flx(i, j, k, to_id + UMX) += fac * flx(i, j, k, from_id + UMX);
          , flx(i, j, k, to_id + UMY) += fac * flx(i, j, k, from_id + UMY);
          , flx(i, j, k, to_id + UMZ) += fac * flx(i, j, k, from_id + UMZ);)
        flx(i, j, k, to_id + UEDEN) += fac * flx(i, j, k, from_id + UEDEN);
        for (int ns = 0; ns < NUM_SPECIES; ++ns) {
          flx(i, j, k, to_id + UFS + ns) += fac * flx(i, j, k, from_id + UFS + ns);
        }
      });
  } // for dir

  amrex::ParallelFor(
    bx, NUM_FIELD, [=] AMREX_GPU_DEVICE(int i, int j, int k, int nfm1) {
      const int nf = 1 + nfm1;
      const Real fac = update_fields ? Real(1.0) / Real(NUM_FIELD) : 1.0;
      const int from_id = update_fields ? nf * NVAR : 0;
      const int to_id = update_fields ? 0 : nf * NVAR;

      dsdt(i, j, k, to_id + URHO) += fac * dsdt(i, j, k, from_id + URHO);
      AMREX_D_TERM(
        dsdt(i, j, k, to_id + UMX) += fac * dsdt(i, j, k, from_id + UMX);
        , dsdt(i, j, k, to_id + UMY) += fac * dsdt(i, j, k, from_id + UMY);
        , dsdt(i, j, k, to_id + UMZ) += fac * dsdt(i, j, k, from_id + UMZ);)
      dsdt(i, j, k, to_id + UEDEN) += fac * dsdt(i, j, k, from_id + UEDEN);
      for (int ns = 0; ns < NUM_SPECIES; ++ns) {
        dsdt(i, j, k, to_id + UFS + ns) += fac * dsdt(i, j, k, from_id + UFS + ns);
      }
    });
}
#endif

  // External source term
  if (do_ext_src) {
    ProbParm const* lprobparm = d_prob_parm;
    const auto geomdata = geom.data();
    const amrex::Real time = state[State_Type].curTime();

    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      fill_ext_src(i, j, k, time, geomdata, sarr, dsdt, *lprobparm);
    });
  }
}
