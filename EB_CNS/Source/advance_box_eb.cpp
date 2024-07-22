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

  const int ncomp = UFA; // UFA because we don't want to change time averages
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
  const bool store_in_vflux = do_diffusion && (NUM_FIELD > 0) && do_vpdf;
  if (store_in_vflux) {
    for (int dir = 0; dir < amrex::SpaceDim; ++dir) {
      vfluxfab[dir].resize(amrex::surroundingNodes(bxg3, dir), NVAR,
                           The_Async_Arena());
      vfluxfab[dir].setVal<RunOn::Device>(0.0);
    }
  }

  // Advance
#if NUM_FIELD > 0
  for (int nf = 1; nf <= NUM_FIELD; ++nf)
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
        amrex::ParallelFor(reconbox, NCHAR,
                           [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) {
                             cns_recon_eb(i, j, k, n, dir, w, wl, wr, recon_scheme,
                                          plm_theta, flag, eb_recon_mode);
                           });
        // 3. Solve Riemann problem for fluxes at cell face
        amrex::ParallelFor(flxbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
          if (!flag(IntVect(AMREX_D_DECL(i, j, k))).isCovered() &&
              !flag(IntVect(AMREX_D_DECL(i, j, k)) -
                    IntVect::TheDimensionVector(dir))
                 .isCovered()) {
            cns_riemann(i, j, k, dir, flx, q, wl, wr, char_sys, recon_char_var);

            bool do_high_order_diff =
              (shock_sensor(i, j, k) < 0.9) &&
              (shock_sensor(IntVect(AMREX_D_DECL(i, j, k)) -
                            IntVect::TheDimensionVector(dir)) < 0.9);
            if (do_high_order_diff) {
              cns_afd_correction_eb(i, j, k, dir, q, flag, flx);
            }
          }
        });
      }

      // Viscous fluxes
      if (do_diffusion) {
        auto const& vflx = store_in_vflux ? vfluxfab[dir].array() : flx;
        amrex::ParallelFor(flxbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
          cns_diff_eb(i, j, k, dir, q, coefs, flag, dxinv, vflx);
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
  for (int nf = 0; nf <= NUM_FIELD; ++nf)
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
    amrex::ParallelFor(bxg3, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
      eb_compute_div(i, j, k, blo, bhi, q, divc, AMREX_D_DECL(fx_in, fy_in, fz_in),
                     AMREX_D_DECL(fx_out, fy_out, fz_out), flag, vfrac, bcent, coefs,
                     AMREX_D_DECL(apx, apy, apz), AMREX_D_DECL(fcx, fcy, fcz), dxinv,
                     do_visc, eb_no_slip, eb_isothermal, eb_wall_temp);
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

  // Zero dsdt for aux variables
  amrex::ParallelFor(bx, NUM_AUX,
                     [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
                       dsdt(i, j, k, UFA + n) = 0.0;
                     });
}