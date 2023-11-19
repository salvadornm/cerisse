#include "CNS.H"
#include "diffusion.H"
#include "hyperbolics.H"
#include "recon.H"
// #include "char_recon.H"

using namespace amrex;

void CNS::compute_dSdt_box(Box const& bx, Array4<const Real>& sarr,
                           Array4<Real>& dsdtarr,
                           std::array<FArrayBox*, AMREX_SPACEDIM> const& flxfab)
{
  BL_PROFILE("CNS::compute_dSdt_box()");

  const Box& bxg2 = amrex::grow(bx, 2);
  const Box& bxg3 = amrex::grow(bx, 3);

  const auto dxinv = geom.InvCellSizeArray();

  // Primitive variables
  FArrayBox qtmp(bxg3, LEN_PRIM, The_Async_Arena());
  Parm const* lparm = d_parm;
  for (int nf = 0; nf <= NUM_FIELD; ++nf) {
    auto const& qfill = qtmp.array(nf * NPRIM);
    amrex::ParallelFor(bxg3, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      cns_ctoprim(i, j, k, nf * NVAR, sarr, qfill, *lparm);
    });
  }

  // Arrays for characteristic reconstruction
  FArrayBox wtmp(bxg3, NCHAR, The_Async_Arena());
  FArrayBox wl_tmp(bxg2, NPRIM, The_Async_Arena());
  FArrayBox wr_tmp(bxg2, NPRIM, The_Async_Arena());
  auto const& w = wtmp.array();
  auto const& wl = wl_tmp.array();
  auto const& wr = wr_tmp.array();

  // Transport coef
  FArrayBox diff_coeff;
  if (do_visc == 1) {
    diff_coeff.resize(bxg2, LEN_COEF, The_Async_Arena());

    for (int nf = 0; nf <= NUM_FIELD; ++nf) {
      auto const& qar_yin = qtmp.const_array(nf * NPRIM + QFS);
      auto const& qar_Tin = qtmp.const_array(nf * NPRIM + QTEMP);
      auto const& qar_rhoin = qtmp.const_array(nf * NPRIM + QRHO);
      auto const& mu = diff_coeff.array(nf * NCOEF + CMU);
      auto const& xi = diff_coeff.array(nf * NCOEF + CXI);
      auto const& lambda = diff_coeff.array(nf * NCOEF + CLAM);
      auto const& rhoD = diff_coeff.array(nf * NCOEF + CRHOD);

      // Get Transport coefs on GPU
      BL_PROFILE("PelePhysics::get_transport_coeffs()");
      auto const* ltransparm = trans_parms.device_trans_parm();
      amrex::launch(bxg2, [=] AMREX_GPU_DEVICE(amrex::Box const& tbx) {
        auto trans = pele::physics::PhysicsType::transport();
        trans.get_transport_coeffs(tbx, qar_yin, qar_Tin, qar_rhoin, rhoD, amrex::Array4<amrex::Real>(), mu, xi,
                                   lambda, ltransparm);
      });

      if (do_les) {
        // Get LES transport coefs
        BL_PROFILE("CNS::LES_transport_coeffs");

        const auto dx = geom.CellSizeArray();
        Real delta = std::pow(AMREX_D_TERM(dx[0], *dx[1], *dx[2]),
                              1.0 / amrex::Real(AMREX_SPACEDIM));

        amrex::ParallelFor(bxg2, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
          Real mu_T, xi_T;
          les_model->mu_T_cc(i, j, k, qar_rhoin, dxinv, delta, Cs, mu_T);
          // les_model->xi_T_cc(i, j, k, qar_rhoin, dxinv, delta, C_I, xi_T);
          xi_T =
            0.0; // TODO: this divu is calculated at cell centre, not face centre

          mu(i, j, k) += mu_T;
          xi(i, j, k) += xi_T;

          for (int n = 0; n < NUM_SPECIES; ++n) rhoD(i, j, k, n) += mu_T / Sc_T;

          Real y[NUM_SPECIES], cp;
          for (int n = 0; n < NUM_SPECIES; ++n) y[n] = qar_yin(i, j, k, n);
          auto eos = pele::physics::PhysicsType::eos();
          eos.RTY2Cp(qar_rhoin(i, j, k), qar_Tin(i, j, k), y, cp);
          lambda(i, j, k) += cp * mu_T / Pr_T;
        });
      }
    }
  }

  // A fab to store the viscous fluxes in VPDF or VSPDF
  FArrayBox vflux_fab;
  // FArrayBox pflux_fab; // unused

  for (int cdir = 0; cdir < AMREX_SPACEDIM; ++cdir) { // Loop through space direction
    const Box& flxbx = amrex::surroundingNodes(bx, cdir);
    const Box& reconbox = amrex::grow(bx, cdir, 1);
    // pflux_fab.resize(flxbx, LEN_STATE);
    // pflux_fab.setVal(0.0);
    if (do_visc == 1) {
      vflux_fab.resize(flxbx, NVAR, The_Async_Arena());
      vflux_fab.setVal(0.0);
    }
    auto const& vflx_arr = vflux_fab.array();

    for (int nf = 0; nf <= NUM_FIELD; ++nf) { // Loop through fields
      // Convert primitive to characteristic
      auto const& q = qtmp.array(nf * NPRIM);
      auto const& flx_arr = flxfab[cdir]->array(nf * NVAR);
      // auto const& p_arr = pflux_fab.array(nf*NVAR);

      // #ifdef CNS_TRUE_CHAR_RECON
      //       // Reconstruction
      //       amrex::ParallelFor(reconbox,
      //       [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
      //         char_recon(i, j, k, cdir, q, wl, wr, recon_scheme, plm_theta,
      //         char_sys, *lparm);
      //       });

      //       // Riemann solver
      //       amrex::ParallelFor(flxbx,
      //       [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
      //         pure_riemann(i, j, k, cdir, flx_arr, wl, wr, *lparm);
      //       });
      // #else
      amrex::ParallelFor(bxg3, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        cns_ctochar(i, j, k, cdir, q, w, char_sys);
      });

      // Reconstruction
      amrex::ParallelFor(
        reconbox, NCHAR, [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
          cns_recon(i, j, k, n, cdir, w, wl, wr, recon_scheme, plm_theta);
        });

      // Solve Riemann problem, store advection and pressure fluxes to flx_arr
      amrex::ParallelFor(flxbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        cns_riemann(i, j, k, cdir, flx_arr, q, wl, wr, char_sys, recon_char_var,
                    *lparm);
      });
      // #endif

      // Store viscous fluxes separately
      if (do_visc) {
        auto const& coefs = diff_coeff.array(nf * NCOEF);
        amrex::ParallelFor(
          flxbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            cns_diff(i, j, k, cdir, q, coefs, dxinv, nf > 0 ? vflx_arr : flx_arr);
          });
      }

    } // for fields

#if NUM_FIELD > 0
    if (do_visc) {
      auto const& flx_arr = flxfab[cdir]->array();
      amrex::ParallelFor(flxbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        // Average the viscous fluxes
        AMREX_D_TERM(vflx_arr(i, j, k, UMX) /= amrex::Real(NUM_FIELD);
                     , vflx_arr(i, j, k, UMY) /= amrex::Real(NUM_FIELD);
                     , vflx_arr(i, j, k, UMZ) /= amrex::Real(NUM_FIELD);)
        vflx_arr(i, j, k, UEDEN) /= amrex::Real(NUM_FIELD);
        for (int n = 0; n < NUM_SPECIES; ++n) {
          vflx_arr(i, j, k, UFS + n) /= amrex::Real(NUM_FIELD);
        }

        amrex::Real mean_rho = sarr(i, j, k, URHO);
        amrex::Real rho;
        // Add to hyperbolic fluxes
        for (int nf = 1; nf <= NUM_FIELD; ++nf) {
          rho = sarr(i, j, k, nf * NVAR + URHO);

          // // <p>*rho/<rho>
          AMREX_D_TERM(flx_arr(i, j, k, nf * NVAR + UMX) +=
                       rho / mean_rho * vflx_arr(i, j, k, UMX);
                       , flx_arr(i, j, k, nf * NVAR + UMY) +=
                         rho / mean_rho * vflx_arr(i, j, k, UMY);
                       , flx_arr(i, j, k, nf * NVAR + UMZ) +=
                         rho / mean_rho * vflx_arr(i, j, k, UMZ);)
          flx_arr(i, j, k, nf * NVAR + UEDEN) +=
            rho / mean_rho * vflx_arr(i, j, k, UEDEN);
          for (int n = 0; n < NUM_SPECIES; ++n) {
            flx_arr(i, j, k, nf * NVAR + UFS + n) +=
              rho / mean_rho * vflx_arr(i, j, k, UFS + n);
          }
        }
      });
    }
#endif
  } // for dir

  // Compute flux divergence
  amrex::ParallelFor(bx, LEN_STATE,
                     [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
                       cns_div(i, j, k, n, dsdtarr,
                               AMREX_D_DECL(flxfab[0]->array(), flxfab[1]->array(),
                                            flxfab[2]->array()),
                               dxinv);
                     });

  // Add external source term
  if (do_ext_src) {
    Parm const* lparm = d_parm;
    ProbParm const* lprobparm = d_prob_parm;
    const auto geomdata = geom.data();
    const amrex::Real time = state[State_Type].curTime();

    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      fill_ext_src(i, j, k, time, geomdata, sarr, dsdtarr, *lparm, *lprobparm);
    });
  }

  Gpu::streamSynchronize();
}
