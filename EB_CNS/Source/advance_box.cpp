#include "CNS.H"
// #include "CNS_hydro_K.H"
#include "CNS_diffusion_K.H"

#include "hyperbolics.H"

using namespace amrex;

/* =============== Pusedocode ===============
 *  if do_visc:       
 *      if vpdf || vspdf || uq:
 *          for nf: 
 *              diff_coef <- get_trans_coef()
 *              flx += diffusion_flux(q)
 *      else:
 *          <diff_coef> <- get_trans_coef()
 *          <flx> += diffusion_flux(<q>)
 * 
 *  if vpdf || vspdf || uq:
 *      for nf: 
 *          q <- ctoprim(S)
 *          w <- qtochar(q)
 *          wL, wR <- recon(w)
 *          flx += pressure_flux(wL, wR)
 *  else:
 *      q <- ctoprim(<S>)
 *      w <- qtochar(q)
 *      wL, wR <- recon(w)  
 *      <flx> += pressure_flux(wL, wR)
 *      <flx> += advection_flux(wL, wR)
 * 
 *  if vspdf || vpdf:
 *      flx <- reduce_average_fields(flx)
 *  else if spdf:
 *      flx <- copy_to_fields(<flx>)
 * 
 *  if vspdf || vpdf || uq:
 *      for nf: 
 *          flx += advection_flux(wL, wR)
 * 
 *  dsdt <- div_and_redistribute(flx)
 *  
 *  if vspdf || vpdf || spdf || uq:
 *      for nf:
 *          dsdt += ext_src(s, t, x)
 *  else:
 *      <dsdt> += ext_src(<s>, t, x)
 *  
 *  if les:
 *      <dsdt> += velocity_les_terms(<s>+<dsdt>*dt)
 *      <dsdt> += species_energy_les_terms(<s>+<dsdt>*dt)
 *  else if spdf:
 *      for nf:
 *          dsdt += velocity_les_terms(<s>+<dsdt>*dt) : this should be the same as s+dsdt*dt
 *  else if vpdf:
 *      for nf:
 *          dsdt += species_energy_les_terms(<s>+<dsdt>*dt)
 * 
 *  if vspdf || vpdf || spdf:
 *      <dsdt> <- reduce_average_fields(dsdt)
 * 
 *  if vspdf || spdf:
 *      for nf: 
 *          dsdt += iem((dsdt-<dsdt>)*dt, dt)
 *  if vspdf || vspdf:
 *      for nf:
 *          dsdt += glm((dsdt-<dsdt>)*dt, dt)
 */

void
CNS::compute_dSdt_box (const Box& bx,
                       Array4<Real const>& sfab,
                       Array4<Real      >& dsdtfab,
                       const std::array<FArrayBox*, AMREX_SPACEDIM>& flxfab
                       //FArrayBox           flxfab[AMREX_SPACEDIM]
                       /*AMREX_D_DECL(
                         Array4<Real    >& fxfab,
                         Array4<Real    >& fyfab,
                         Array4<Real    >& fzfab)*/)
{
    BL_PROFILE("CNS::compute_dSdt_box()");
    
    const Box& bxg1 = amrex::grow(bx,1);
    const Box& bxg2 = amrex::grow(bx,2);
    const Box& bxg3 = amrex::grow(bx,3);
    
    const auto dxinv = geom.InvCellSizeArray();
    
    // Primitive variables
    FArrayBox qtmp(bxg3, LEN_PRIM);
    auto const& q = qtmp.array();

    Parm const* lparm = d_parm;
    // for (int nf = 0; nf <= NUM_FIELD; ++nf){
    //     auto const& qfill = qtmp.array(nf*NPRIM);
    amrex::ParallelFor(bxg3,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        cns_ctoprim(i, j, k, sfab, q, *lparm);
    });
    // }

    // Slopes
    // FArrayBox slopetmp(bxg1, 5+NUM_SPECIES);
    // auto const& slope = slopetmp.array();
    FArrayBox wtmp(bxg3, NCHAR);
    FArrayBox wl_tmp(bxg2, NCHAR);
    FArrayBox wr_tmp(bxg2, NCHAR);
    auto const& w = wtmp.array();
    auto const& wl = wl_tmp.array();
    auto const& wr = wr_tmp.array();

    // Transport coef
    FArrayBox diff_coeff;
    if (do_visc == 1) {
        diff_coeff.resize(bxg2, NCOEF);
       
        auto const& qar_yin   = qtmp.array(QFS);
        auto const& qar_Tin   = qtmp.array(QTEMP);
        auto const& qar_rhoin = qtmp.array(QRHO);
        auto const& rhoD   = diff_coeff.array(CRHOD);
        auto const& mu     = diff_coeff.array(CMU);
        auto const& xi     = diff_coeff.array(CXI);
        auto const& lambda = diff_coeff.array(CLAM);
        
        BL_PROFILE("PelePhysics::get_transport_coeffs()");
        // Get Transport coefs on GPU.
        auto const* ltransparm = trans_parms.device_trans_parm();
        amrex::launch(bxg2, [=] AMREX_GPU_DEVICE(amrex::Box const& tbx) {
            auto trans = pele::physics::PhysicsType::transport();
            trans.get_transport_coeffs(
                tbx, qar_yin, qar_Tin, qar_rhoin, rhoD, mu, xi, lambda, ltransparm);
        });
    }

    // Loop directions
    for (int cdir = 0; cdir < AMREX_SPACEDIM; ++cdir) {
        // Convert primitive to characteristic
        amrex::ParallelFor(bxg3,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            cns_ctochar(i, j, k, cdir, q, w);
        });

        // Reconstruction
        const Box& qlrbox = amrex::grow(bx,cdir,1);
        amrex::ParallelFor(qlrbox, NCHAR,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept {
            cns_recon(i, j, k, n, cdir, w, wl, wr, recon_scheme, plm_theta);
        });

        // Solve Riemann problem
        const Box& flxbx = amrex::surroundingNodes(bx,cdir);
        auto const& flx_arr = flxfab[cdir]->array();
        amrex::ParallelFor(flxbx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            cns_riemann(i, j, k, cdir, flx_arr, q, wl, wr, *lparm);
            for (int n = UFA; n < NUM_AUX; ++n) flx_arr(i,j,k,n) = Real(0.0);
        });

        // Get viscous fluxes
        if (do_visc == 1) {
            auto const& coefs = diff_coeff.array();
            amrex::ParallelFor(flxbx,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                cns_diff(i, j, k, cdir, q, coefs, dxinv, flx_arr);
            });
        }
    }
    
    amrex::ParallelFor(bx, NVAR,
    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        cns_div(i, j, k, n, dsdtfab, AMREX_D_DECL(flxfab[0]->array(),flxfab[1]->array(),flxfab[2]->array()), dxinv);
    });

//     if (do_ext_src) {
//     ...
//         const Real g = gravity;
//         const int irho = URHO;
// #if (AMREX_SPACEDIM == 2)
//         const int imz = UMY;
// #elif (AMREX_SPACEDIM == 3)
//         const int imz = UMZ;
// #endif
//         const int irhoE = UEDEN;
//         amrex::ParallelFor(bx,
//         [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
//         {
//             dsdtfab(i,j,k,imz  ) += g * sfab(i,j,k,irho);
//             dsdtfab(i,j,k,irhoE) += g * sfab(i,j,k,imz);
//         });
//     }

    Gpu::streamSynchronize();
}
