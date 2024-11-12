#include "CNS.H"
#include "central_scheme.H"
#include "diffusion_eb.H"
#include "hydro.H"
#include "hydro_eb.H"
#include "recon_eb.H"

#if (AMREX_SPACEDIM == 2)
#include <AMReX_EBMultiFabUtil_2D_C.H>
#elif (AMREX_SPACEDIM == 3)
#include <AMReX_EBMultiFabUtil_3D_C.H>
#endif
#include <AMReX_EBFArrayBox.H>
#include <AMReX_EB_Redistribution.H>
#include <AMReX_MultiCutFab.H>

using namespace amrex;

void CNS::compute_dSdt_box_eb(
  const Box& bx, Array4<const Real> const& sarr, Array4<Real> const& dsdt,
  std::array<FArrayBox*, AMREX_SPACEDIM> const& flxfab,
  Array4<const EBCellFlag> const& flag, Array4<const Real> const& vfrac,
  AMREX_D_DECL(Array4<const Real> const& apx, Array4<const Real> const& apy,
               Array4<const Real> const& apz),
  AMREX_D_DECL(Array4<const Real> const& fcx, Array4<const Real> const& fcy,
               Array4<const Real> const& fcz),
  Array4<const Real> const& bcent, int as_crse, Array4<Real> const& rr_drho_crse,
  Array4<const int> const& rr_flag_crse, int as_fine, Array4<Real> const& dm_as_fine,
  Array4<const int> const& lev_mask, Real dt, Array4<const Real>& shock_sensor)
{
  BL_PROFILE("CNS::compute_dSdt_box_eb()");
  
#if NUM_FIELD > 0
  // auto const& geomdata = geom.data();
  // const Real x = geomdata.ProbLo(0) + (bx.bigEnd(0) + 0.5) * geomdata.CellSize(0);
  // const bool update_fields = (x > -1.0) && (x < 5.0); // Do SF or just the mean? TODO: make this a tagging
  const bool update_fields = true;
#else
  const bool update_fields = false;
#endif
  const int ncomp = update_fields ? UFA : NVAR;
  const Box& bxg3 = amrex::grow(bx, 3);
  const Box& bxg4 = amrex::grow(bx, 4);
  const Box& bxg5 = amrex::grow(bx, 5);
  const Box& bxg6 = amrex::grow(bx, 6);
  const auto dx = geom.CellSizeArray();
  const auto dxinv = geom.InvCellSizeArray();
  const bool do_diffusion = do_visc || do_les || buffer_box.ok();

  // Prepare FABs to store data
  FArrayBox divcfab(bxg3, ncomp, The_Async_Arena()); // For redistribution
  FArrayBox qfab(bxg6, LEN_PRIM, The_Async_Arena()); // Primitive variables
  FArrayBox wfab(bxg6, NCHAR, The_Async_Arena());    // FABs for reconstruction
  FArrayBox wlfab(bxg5, NCHAR, The_Async_Arena());
  FArrayBox wrfab(bxg5, NCHAR, The_Async_Arena());
  FArrayBox coefsfab; // Diffusion coefficients
  if (do_diffusion) { coefsfab.resize(bxg4, LEN_COEF, The_Async_Arena()); }
  auto const& w = wfab.array();
  auto const& wl = wlfab.array();
  auto const& wr = wrfab.array();

  // Temporary flux before redistribution
  std::array<FArrayBox, amrex::SpaceDim> flux_tmp;
  for (int dir = 0; dir < amrex::SpaceDim; ++dir) {
    flux_tmp[dir].resize(amrex::surroundingNodes(bxg3, dir), ncomp,
                         The_Async_Arena());
    flux_tmp[dir].setVal<RunOn::Device>(0.);
  }
  // Store viscous fluxes separately in V/VSPDF
  std::array<FArrayBox, amrex::SpaceDim> vfluxfab;

  const bool store_in_vflux = do_diffusion && do_vpdf && update_fields;  
  if (store_in_vflux) {
    for (int dir = 0; dir < amrex::SpaceDim; ++dir) {
      vfluxfab[dir].resize(amrex::surroundingNodes(bxg3, dir), NVAR,
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
    auto const& q = qfab.array(nf * NPRIM);
    auto const& coefs = coefsfab.array(nf * NCOEF);

    // Prepare primitive variables
    amrex::ParallelFor(bxg6, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
      if (!flag(i, j, k).isCovered()) { cns_ctoprim(i, j, k, nf * NVAR, sarr, q); }
    });

    // Prepare transport coefs
    if (do_diffusion) {
      auto const& qar_yin = qfab.const_array(nf * NPRIM + QFS);
      auto const& qar_Tin = qfab.const_array(nf * NPRIM + QTEMP);
      auto const& qar_rhoin = qfab.const_array(nf * NPRIM + QRHO);
      auto const& mu = coefsfab.array(nf * NCOEF + CMU);      // dynamic viscosity
      auto const& xi = coefsfab.array(nf * NCOEF + CXI);      // bulk viscosity
      auto const& lambda = coefsfab.array(nf * NCOEF + CLAM); // thermal conductivity
      auto const& rhoD = coefsfab.array(nf * NCOEF + CRHOD);  // rho * diffusivity

      // Physical transport coefs
      if (do_visc) {
        BL_PROFILE("PelePhysics::get_transport_coeffs()");
        auto const* ltransparm = trans_parms.device_trans_parm();
        amrex::ParallelFor(bxg4, [=](int i, int j, int k) {
          if (!flag(i, j, k).isCovered()) {
            Real muloc, xiloc, lamloc, Ddiag[NUM_SPECIES];
            Real* dummy_chi_mix = nullptr;
            Real yin[NUM_SPECIES];
            for (int n = 0; n < NUM_SPECIES; ++n) { yin[n] = qar_yin(i, j, k, n); }

            auto trans = pele::physics::PhysicsType::transport();
            const bool wtr_get_xi = true;
            const bool wtr_get_mu = true;
            const bool wtr_get_lam = true;
            const bool wtr_get_Ddiag = true;
            const bool wtr_get_chi = false;
            trans.transport(wtr_get_xi, wtr_get_mu, wtr_get_lam, wtr_get_Ddiag,
                            wtr_get_chi, qar_Tin(i, j, k), qar_rhoin(i, j, k), yin,
                            Ddiag, dummy_chi_mix, muloc, xiloc, lamloc, ltransparm);

            // Copy to output
            mu(i, j, k) = muloc;
            xi(i, j, k) = xiloc;
            lambda(i, j, k) = lamloc;
            for (int n = 0; n < NUM_SPECIES; ++n) { rhoD(i, j, k, n) = Ddiag[n]; }
          }
        });
      }

      // LES diffucsion coefs
      if (do_les) { amrex::Abort("LES for EB not implemented yet"); }

      // Buffer region
      if (buffer_box.ok()) {
        const auto problo = geom.ProbLo();
        amrex::ParallelFor(bxg4, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
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
    }

    // Compute fluxes for each space direction
    for (int dir = 0; dir < amrex::SpaceDim; ++dir) {
      const Box& flxbx = amrex::surroundingNodes(bxg3, dir);
      auto const& flx = flux_tmp[dir].array(nf * NVAR);

      // Hydro/hyperbolic fluxes
      if (do_hydro) {
        // Shock-capturing scheme
        // 1. Convert primitive to characteristic at cell centre
        const Box& charbox = amrex::grow(bxg3, dir, 3);
        amrex::ParallelFor(charbox, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
          if (!flag(i, j, k).isCovered()) {
            cns_ctochar(i, j, k, dir, q, w, char_sys);
          }
        });
        // 2. FD interpolation to cell face
        const Box& reconbox = amrex::grow(bxg3, dir, 1);
        amrex::ParallelFor(
          TypeList<CompileTimeOptions<1, 2, 3, 4, 5, 6>, CompileTimeOptions<0, 1>>{},
          {recon_scheme, eb_recon_mode}, reconbox, NCHAR,
          [=] AMREX_GPU_DEVICE(int i, int j, int k, int n,
                               auto captured_recon_scheme,
                               auto captured_eb_recon_mode) {
            cns_recon_eb<captured_recon_scheme, captured_eb_recon_mode>(
              i, j, k, n, dir, w, wl, wr, plm_theta, flag);
          });
        // 3. Solve Riemann problem for fluxes at cell face
        amrex::ParallelFor(flxbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
          if (!flag(IntVect(AMREX_D_DECL(i, j, k))).isCovered() &&
              !flag(IntVect(AMREX_D_DECL(i, j, k)) -
                    IntVect::TheDimensionVector(dir))
                 .isCovered()) {
            cns_riemann(i, j, k, dir, flx, q, wl, wr, char_sys, recon_char_var);

            // bool do_high_order_diff =
            //   (shock_sensor(i, j, k) < 0.95) &&
            //   (shock_sensor(IntVect(AMREX_D_DECL(i, j, k)) -
            //                 IntVect::TheDimensionVector(dir)) < 0.95);
            // if (do_high_order_diff) {
            //   cns_afd_correction_eb(i, j, k, dir, q, flag, flx);
            // }
          }
        });
      }

      // Viscous fluxes
      if (do_diffusion) {
        auto const& vflx = store_in_vflux ? vfluxfab[dir].array() : flx;
        amrex::ParallelFor(flxbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
          const IntVect iv{AMREX_D_DECL(i, j, k)};
          const IntVect ivm(iv - IntVect::TheDimensionVector(dir));          
          if (flag(iv).isCovered() || flag(ivm).isCovered()) { return; }

          Real flx_tmp[NVAR] = {0.0};
          cns_diff_eb(iv, dir, q, coefs, flag, dxinv, flx_tmp);

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
                  // flx_tmp[UEDEN] = qw * lohi; // TODO: no wall model heat flux for now
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
      const Box& flxbx = amrex::surroundingNodes(bxg3, dir);
      auto const& flx = flux_tmp[dir].array();
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

  ///////////////// PART DIFFERENT FROM NON-EB VERSION /////////////////

  // Compute flux divergence for EB
#if NUM_FIELD > 0
  for (int nf = nf_start; nf <= nf_end; ++nf)
  // #else const int nf = 0; declared before
#endif
  {
    // These fluxes are computed above, they are on face centres
    AMREX_D_TERM(auto const& fx_in = flux_tmp[0].array(nf * NVAR);
                 , auto const& fy_in = flux_tmp[1].array(nf * NVAR);
                 , auto const& fz_in = flux_tmp[2].array(nf * NVAR);)

    // These fluxes are on face centroids, they are calculated in eb_compute_div
    // and are the fluxes that go into the flux registers
    AMREX_D_TERM(auto const& fx_out = flxfab[0]->array(nf * NVAR);
                 , auto const& fy_out = flxfab[1]->array(nf * NVAR);
                 , auto const& fz_out = flxfab[2]->array(nf * NVAR);)

    // Because we will do redistribute, we put the divergence into divc rather than
    // directly into dsdt
    auto const& divc = divcfab.array(nf * NVAR);
    auto const& q = qfab.array(nf * NPRIM);
    auto const& coefs = coefsfab.array(nf * NCOEF);

    auto const& blo = bx.smallEnd();
    auto const& bhi = bx.bigEnd();

    // This does the divergence and cut face wall BCs
    const auto problo = geom.ProbLo();
    amrex::ParallelFor(bxg3, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
      const Real x = (i + 0.5) * dx[0] + problo[0];
      eb_compute_div(i, j, k, blo, bhi, q, divc, AMREX_D_DECL(fx_in, fy_in, fz_in),
                     AMREX_D_DECL(fx_out, fy_out, fz_out), flag, vfrac, bcent, coefs,
                     AMREX_D_DECL(apx, apy, apz), AMREX_D_DECL(fcx, fcy, fcz), dxinv,
                     do_hydro, do_visc, eb_no_slip, eb_isothermal, eb_wall_temp,
                     eb_wall_model && (x < -1.0));
    });
  } // for fields

  // Now do redistribution
  {
    BL_PROFILE("Redistribution");
    FArrayBox redistwgt_fab(bxg3, 1);        // redistribution weights
    FArrayBox srd_update_scale_fab(bxg3, 1); // redistribution weights
    auto const& redistwgt = redistwgt_fab.array();
    auto const& srd_update_scale = srd_update_scale_fab.array();
    amrex::ParallelFor(bxg3, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      redistwgt(i, j, k) = vfrac(i, j, k);
      srd_update_scale(i, j, k) =
        eb_weight; // = 1.0 more stable for detondiff, = 0.5 for airfoil. Why?
    });

    bool use_wts_in_divnc = true;
    int srd_max_order = 2;
    Real target_volfrac = 0.5;

    FArrayBox tmpfab(bxg3, ncomp, The_Async_Arena());
    Array4<Real> scratch = tmpfab.array();
    if (redistribution_type == "FluxRedist") {
      amrex::ParallelFor(bxg3, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        scratch(i, j, k) = redistwgt(i, j, k);
      });
    }

    auto const& divc = divcfab.array();
    int level_mask_not_covered = CNSConstants::level_mask_notcovered;
    Real fac_for_deltaR = 1.0;
    amrex::ApplyMLRedistribution(
      bx, ncomp, dsdt, divc, sarr, scratch, flag, AMREX_D_DECL(apx, apy, apz), vfrac,
      AMREX_D_DECL(fcx, fcy, fcz), bcent, &phys_bc, geom, dt, redistribution_type,
      as_crse, rr_drho_crse, rr_flag_crse, as_fine, dm_as_fine, lev_mask,
      level_mask_not_covered, fac_for_deltaR, use_wts_in_divnc, 0, srd_max_order,
      target_volfrac, srd_update_scale);
  }

  //////////////////////// END EB CALCULATIONS ////////////////////////

#if NUM_FIELD > 0
{
  // BL_PROFILE("New_field_averaging_EB");
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

//   // Copy mean field to all fields if not doing SF
//   if (!update_fields) {
//     for (int dir = 0; dir < amrex::SpaceDim; ++dir) {
//       const Box& flxbx = amrex::surroundingNodes(bx, dir);
//       auto const& flx = flxfab[dir]->array();

//       amrex::ParallelFor(
//         flxbx, NUM_FIELD, [=] AMREX_GPU_DEVICE(int i, int j, int k, int nfm1) {
//           const int nf = 1 + nfm1;
//           flx(i, j, k, nf * NVAR + URHO) = flx(i, j, k, URHO);
//           AMREX_D_TERM(flx(i, j, k, nf * NVAR + UMX) = flx(i, j, k, UMX);
//                        , flx(i, j, k, nf * NVAR + UMY) = flx(i, j, k, UMY);
//                        , flx(i, j, k, nf * NVAR + UMZ) = flx(i, j, k, UMZ);)
//           flx(i, j, k, nf * NVAR + UEDEN) = flx(i, j, k, UEDEN);
//           for (int ns = 0; ns < NUM_SPECIES; ++ns) {
//             flx(i, j, k, nf * NVAR + UFS + ns) = flx(i, j, k, UFS + ns);
//           }
//         });
//     } // for dir

//     amrex::ParallelFor(
//       bx, NUM_FIELD, [=] AMREX_GPU_DEVICE(int i, int j, int k, int nfm1) {
//         const int nf = 1 + nfm1;
//         dsdt(i, j, k, nf * NVAR + URHO) = dsdt(i, j, k, URHO);
//         AMREX_D_TERM(dsdt(i, j, k, nf * NVAR + UMX) = dsdt(i, j, k, UMX);
//                      , dsdt(i, j, k, nf * NVAR + UMY) = dsdt(i, j, k, UMY);
//                      , dsdt(i, j, k, nf * NVAR + UMZ) = dsdt(i, j, k, UMZ);)
//         dsdt(i, j, k, nf * NVAR + UEDEN) = dsdt(i, j, k, UEDEN);
//         for (int ns = 0; ns < NUM_SPECIES; ++ns) {
//           dsdt(i, j, k, nf * NVAR + UFS + ns) = dsdt(i, j, k, UFS + ns);
//         }
//       });
//   } else {
//     // Average to mean
//     for (int dir = 0; dir < amrex::SpaceDim; ++dir) {
//       const Box& flxbx = amrex::surroundingNodes(bx, dir);
//       auto const& flx = flxfab[dir]->array();

//       amrex::ParallelFor(
//         flxbx, NUM_FIELD, [=] AMREX_GPU_DEVICE(int i, int j, int k, int nfm1) {
//           const int nf = 1 + nfm1;
//           const Real invNF = Real(1.0) / Real(NUM_FIELD);
//           flx(i, j, k, URHO) += invNF * flx(i, j, k, nf * NVAR + URHO);
//           AMREX_D_TERM(flx(i, j, k, UMX) += invNF * flx(i, j, k, nf * NVAR + UMX);
//                        , flx(i, j, k, UMY) += invNF * flx(i, j, k, nf * NVAR + UMY);
//                        , flx(i, j, k, UMZ) += invNF * flx(i, j, k, nf * NVAR + UMZ);)
//           flx(i, j, k, UEDEN) += invNF * flx(i, j, k, nf * NVAR + UEDEN);
//           for (int ns = 0; ns < NUM_SPECIES; ++ns) {
//             flx(i, j, k, UFS + ns) += invNF * flx(i, j, k, nf * NVAR + UFS + ns);
//           }          
//         });
//     } // for dir
//   }
}
#endif

  // External source term
  if (do_ext_src) {
    ProbParm const* lprobparm = d_prob_parm;
    const auto geomdata = geom.data();
    const amrex::Real time = state[State_Type].curTime();

    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      if (!flag(i, j, k).isCovered()) {
        fill_ext_src(i, j, k, time, geomdata, sarr, dsdt, *lprobparm);
      }
    });
  }
}
