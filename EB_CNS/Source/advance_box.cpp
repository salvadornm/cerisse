#include "CNS.H"
#include "central_scheme.H"
#include "diffusion.H"
#include "hydro.H"
#include "recon.H"

using namespace amrex;

void CNS::compute_dSdt_box(Box const& bx, Array4<const Real>& sarr,
                           Array4<Real>& dsdt,
                           std::array<FArrayBox*, AMREX_SPACEDIM> const& flxfab,
                           Array4<const Real>& shock_sensor)
{
  BL_PROFILE("CNS::compute_dSdt_box()");

  const int ncomp = UFA; // UFA because we don't want to change time averages
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
  const bool store_in_vflux = do_diffusion && (NUM_FIELD > 0) && do_vpdf;
  if (store_in_vflux) {
    for (int dir = 0; dir < amrex::SpaceDim; ++dir) {
      vfluxfab[dir].resize(amrex::surroundingNodes(bx, dir), NVAR,
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
        amrex::Array4<amrex::Real> chi; // dummy Soret effect coef
        auto const* ltransparm = trans_parms.device_trans_parm();
        amrex::launch(bxg2, [=] AMREX_GPU_DEVICE(amrex::Box const& tbx) {
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
          amrex::ParallelFor(
            reconbox, NCHAR, [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) {
              cns_recon(i, j, k, n, dir, w, wl, wr, recon_scheme, plm_theta);
            });
          // 3. Solve Riemann problem for fluxes at cell face
          amrex::ParallelFor(flxbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            cns_riemann(i, j, k, dir, flx, q, wl, wr, char_sys, recon_char_var);

            bool do_high_order_diff =
              (shock_sensor(i, j, k) < 0.9) &&
              (shock_sensor(IntVect(AMREX_D_DECL(i, j, k)) -
                            IntVect::TheDimensionVector(dir)) < 0.9);
            if (do_high_order_diff) { cns_afd_correction(i, j, k, dir, q, flx); }
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
          cns_diff(i, j, k, dir, q, coefs, dxinv, vflx);
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