#include "CNS.H"

#include <AMReX_EBFArrayBox.H>
#include <AMReX_MultiCutFab.H>
#if (AMREX_SPACEDIM == 2)
#include <AMReX_EBMultiFabUtil_2D_C.H>
#elif (AMREX_SPACEDIM == 3)
#include <AMReX_EBMultiFabUtil_3D_C.H>
#endif

#include "hyperbolics.H"
#include "recon_eb.H"
#include "CNS_divop_K.H"
#include "CNS_diffusion_eb_K.H"

using namespace amrex;

void
CNS::compute_dSdt_box_eb (const Box& bx,
                          Array4<const Real> const& sarr,
                          Array4<      Real> const& dsdtarr,
                          std::array<FArrayBox*, AMREX_SPACEDIM> const& flxfab,
                          Array4<const EBCellFlag> const& flag,
                          Array4<const Real> const& vfrac,
                          AMREX_D_DECL(
                            Array4<const Real> const& apx,
                            Array4<const Real> const& apy,
                            Array4<const Real> const& apz),
                          AMREX_D_DECL(
                            Array4<const Real> const& fcx,
                            Array4<const Real> const& fcy,
                            Array4<const Real> const& fcz),
                          Array4<const Real> const& bcent,
                          int as_crse,
                          Array4<Real> const& drho_as_crse,
                          Array4<const int> const& rrflag_as_crse,
                          int as_fine,
                          Array4<Real> const& dm_as_fine,
                          Array4<const int> const& lev_mask,
                          Real dt)
{
  BL_PROFILE("CNS::compute_dSdt_box_eb()");

  const Box& bxg1 = amrex::grow(bx,1);
  const Box& bxg2 = amrex::grow(bx,2);
  const Box& bxg3 = amrex::grow(bx,3);
  const Box& bxg4 = amrex::grow(bx,4);
  const Box& bxg5 = amrex::grow(bx,5);
  const Box& bxg6 = amrex::grow(bx,6);

  const auto dxinv = geom.InvCellSizeArray();

  // Quantities for redistribution
  FArrayBox divc, redistwgt;
  divc.resize(bxg2, LEN_STATE); //div before redist
  divc.setVal<RunOn::Device>(0.0);

  // Primitive variables
  FArrayBox qtmp(bxg6, LEN_PRIM);
  qtmp.setVal<RunOn::Device>(0.0); //set covered cells to zero to prevent error in compute_diff_wallflux
  Parm const* lparm = d_parm;
  for (int nf = 0; nf <= NUM_FIELD; ++nf){
    auto const& qfill = qtmp.array(nf*NPRIM);
    amrex::ParallelFor(bxg6,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
      if (!flag(i,j,k).isCovered())
        cns_ctoprim(i, j, k, nf*NVAR, sarr, qfill, *lparm);
    });
  }

  // Arrays for characteristic reconstruction
  FArrayBox wtmp(bxg6, NCHAR);
  FArrayBox wl_tmp(bxg5, NPRIM);
  FArrayBox wr_tmp(bxg5, NPRIM);
  auto const& w = wtmp.array();
  auto const& wl = wl_tmp.array();
  auto const& wr = wr_tmp.array();

  // Transport coef
  FArrayBox diff_coeff;
  if (do_visc == 1) {
    diff_coeff.resize(bxg4, LEN_COEF);

    for (int nf = 0; nf <= NUM_FIELD; ++nf){
      auto const& qar_yin   = qtmp.array(nf*NPRIM + QFS);
      auto const& qar_Tin   = qtmp.array(nf*NPRIM + QTEMP);
      auto const& qar_rhoin = qtmp.array(nf*NPRIM + QRHO);
      auto const& mu     = diff_coeff.array(nf*NCOEF + CMU);
      auto const& xi     = diff_coeff.array(nf*NCOEF + CXI);
      auto const& lambda = diff_coeff.array(nf*NCOEF + CLAM);
      auto const& rhoD   = diff_coeff.array(nf*NCOEF + CRHOD);

      // Get Transport coefs on GPU
      BL_PROFILE("PelePhysics::get_transport_coeffs()");      
      auto const* ltransparm = trans_parms.device_trans_parm();
      amrex::launch(bxg4, [=] AMREX_GPU_DEVICE(amrex::Box const& tbx) {
        auto trans = pele::physics::PhysicsType::transport();
        trans.get_transport_coeffs(
            tbx, qar_yin, qar_Tin, qar_rhoin, rhoD, mu, xi, lambda, ltransparm);
      });
    }
  }

  // Initialize dm_as_fine to 0
  if (as_fine) {
    amrex::ParallelFor(bxg1, LEN_STATE,
    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept {
      dm_as_fine(i,j,k,n) = 0.;
    });
  }

  // A fab to store the viscous fluxes in VPDF or VSPDF
  FArrayBox vflux_fab, pflux_fab;

  // EB weights
  GpuArray<const Real, 3> weights{0.0, 1.0, 0.5};

  // Temporary flux before redistribution
  FArrayBox flux_tmp[AMREX_SPACEDIM];

  for (int cdir = 0; cdir < AMREX_SPACEDIM; ++cdir) { // Loop through space direction
    const Box& flxbx = amrex::surroundingNodes(bxg3,cdir);
    const Box& reconbox = amrex::grow(bxg3,cdir,1);
    // Tmp fab for viscous fluxes
    if (do_visc == 1) {
      vflux_fab.resize(flxbx, NVAR);
      vflux_fab.setVal(0.0);            
    }
    auto const& vflx_arr = vflux_fab.array();
    // Tmp fab for pressure fluxes (not used)
    pflux_fab.resize(flxbx, NVAR);
    pflux_fab.setVal(0.0);
    auto const& p_arr = pflux_fab.array();
    // Tmp fab for all fluxes
    flux_tmp[cdir].resize(flxbx,LEN_STATE);
    flux_tmp[cdir].setVal<RunOn::Device>(0.);

    for (int nf = 0; nf <= NUM_FIELD; ++nf) { // Loop through fields
      // Convert primitive to characteristic
      auto const& q = qtmp.array(nf*NPRIM);
      amrex::ParallelFor(bxg6,
      [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        if (!flag(i,j,k).isCovered())
          cns_ctochar(i, j, k, cdir, q, w, 0);
      });

      // Reconstruction            
      amrex::ParallelFor(reconbox, NCHAR,
      [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept {
        cns_recon_eb(i, j, k, n, cdir, w, wl, wr, recon_scheme, plm_theta, flag, 1); // Type 1 EB
      });

      // Solve Riemann problem, store advection and pressure fluxes to flx_arr
      auto const& flx_arr = flux_tmp[cdir].array(nf*NVAR);      
      amrex::ParallelFor(flxbx,
      [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        if ( !flag(amrex::IntVect(AMREX_D_DECL(i,j,k))).isCovered() && 
             !flag(amrex::IntVect(AMREX_D_DECL(i,j,k))-amrex::IntVect::TheDimensionVector(cdir)).isCovered() )
          cns_riemann(i, j, k, cdir, flx_arr, p_arr, q, wl, wr, 0, *lparm);
      });

      // Store viscous fluxes separately
      if (do_visc == 1) {
        auto const& coefs = diff_coeff.array(nf*NCOEF);
        amrex::ParallelFor(flxbx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
          cns_diff_eb(i, j, k, cdir, q, coefs, flag, dxinv, weights, nf > 0 ? vflx_arr : flx_arr);
        });
      }
      
    } //loop for fields

#if (NUM_FIELD > 0)
    if (do_visc == 1) {
      auto const& flx_arr = flxfab[cdir]->array();
      amrex::ParallelFor(flxbx,
      [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        // Average the viscous fluxes
        AMREX_D_TERM(
            vflx_arr(i, j, k, UMX) /= amrex::Real(NUM_FIELD); ,
            vflx_arr(i, j, k, UMY) /= amrex::Real(NUM_FIELD); ,
            vflx_arr(i, j, k, UMZ) /= amrex::Real(NUM_FIELD);)
        vflx_arr(i, j, k, UEDEN) /= amrex::Real(NUM_FIELD);
        for (int n = 0; n < NUM_SPECIES; ++n) {
          vflx_arr(i, j, k, UFS+n) /= amrex::Real(NUM_FIELD);
        }

        // Add to hyperbolic fluxes of all fields
        amrex::Real mean_rho = sarr(i, j, k, URHO);
        amrex::Real rho;
        for (int nf = 1; nf <= NUM_FIELD; ++nf) {
          rho = sarr(i, j, k, nf*NVAR + URHO);
          
          AMREX_D_TERM(
              flx_arr(i, j, k, nf*NVAR + UMX) += rho / mean_rho * vflx_arr(i, j, k, UMX); ,
              flx_arr(i, j, k, nf*NVAR + UMY) += rho / mean_rho * vflx_arr(i, j, k, UMY); ,
              flx_arr(i, j, k, nf*NVAR + UMZ) += rho / mean_rho * vflx_arr(i, j, k, UMZ);)            
          flx_arr(i, j, k, nf*NVAR + UEDEN) += rho / mean_rho * vflx_arr(i, j, k, UEDEN);      
          for (int n = 0; n < NUM_SPECIES; ++n) {
            flx_arr(i, j, k, nf*NVAR + UFS+n) += rho / mean_rho * vflx_arr(i, j, k, UFS+n);
          }
        }
      });
    }
#endif
  } //loop for dir

  ///////////////// PART DIFFERENT FROM NON-EB VERSION /////////////////

  auto const& blo = bx.smallEnd();
  auto const& bhi = bx.bigEnd();

  // Compute flux divergence
  for (int nf = 0; nf <= NUM_FIELD; ++nf) { // Loop through fields
    // These are the fluxes we computed above -- they live at face centers
    AMREX_D_TERM(auto const& fx_in_arr = flux_tmp[0].array(nf*NVAR);,
                 auto const& fy_in_arr = flux_tmp[1].array(nf*NVAR);,
                 auto const& fz_in_arr = flux_tmp[2].array(nf*NVAR););

    // These are the fluxes on face centroids -- they are defined in eb_compute_div
    // and are the fluxes that go into the flux registers
    AMREX_D_TERM(auto const& fx_out_arr = flxfab[0]->array(nf*NVAR);,
                 auto const& fy_out_arr = flxfab[1]->array(nf*NVAR);,
                 auto const& fz_out_arr = flxfab[2]->array(nf*NVAR););

    // Because we will do redistribute, we put the divergence into divc
    // rather than directly into dsdtfab
    auto const& divc_arr = divc.array(nf*NVAR);

    auto const& q = qtmp.array(nf*NPRIM);
    auto const& coefs = diff_coeff.array(nf*NCOEF);

    amrex::ParallelFor(bxg2,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
      // This does the divergence but not the redistribution -- we will do that later
      // We do compute the weights here though (WHY??? WHY NOT COMPUTE LATER? NO LOOP THERE)
      eb_compute_div(i, j, k, blo, bhi, q, divc_arr,
                     AMREX_D_DECL(fx_in_arr, fy_in_arr, fz_in_arr),
                     AMREX_D_DECL(fx_out_arr, fy_out_arr, fz_out_arr),
                     flag, vfrac, bcent, coefs,
                     AMREX_D_DECL(apx, apy, apz),
                     AMREX_D_DECL(fcx, fcy, fcz), dxinv, *lparm, do_visc, 
                     eb_no_slip, eb_isothermal, eb_wall_temp);
    });
  
    // Now do redistribution
    cns_flux_redistribute(bx, dsdtarr, divc_arr, q, vfrac, flag,
                          as_crse, drho_as_crse, rrflag_as_crse, 
                          as_fine, dm_as_fine, lev_mask, dt);
  } //loop for fields

  //////////////////////// END EB CALCULATIONS ////////////////////////

  // Add external source term
  if (do_ext_src) {
    Parm const* lparm = d_parm;
    ProbParm const* lprobparm = d_prob_parm;
    const auto geomdata = geom.data();
    const amrex::Real time = state[State_Type].curTime();

    amrex::ParallelFor(bx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
      fill_ext_src(i, j, k, time, geomdata, sarr, dsdtarr, *lparm, *lprobparm);
    });
  }

  Gpu::streamSynchronize();
}
