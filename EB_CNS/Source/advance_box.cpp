#include <CNS.H>
#include <CNS_hydro_K.H>
#include <CNS_diffusion_K.H>
#include "diffcoef.H"

#include <AMReX_EBFArrayBox.H>
#include <AMReX_MultiCutFab.H>

using namespace amrex;

void
CNS::compute_dSdt_box(const Box& bx,
                      Array4<Real const>& sfab,
                      Array4<Real      >& dsdtfab,
                      /*const std::array<FArrayBox*, AMREX_SPACEDIM>& flux*/
                      AMREX_D_DECL(
                        Array4<Real    >& fxfab,
                        Array4<Real    >& fyfab,
                        Array4<Real    >& fzfab))
{
    BL_PROFILE("CNS::compute_dSdt__box()");

    const auto dxinv = geom.InvCellSizeArray();

    FArrayBox qtmp, slopetmp, diff_coeff;

    Parm const* lparm = d_parm;

    // AMREX_D_TERM(auto const& fxfab = flux[0]->array();,
    //              auto const& fyfab = flux[1]->array();,
    //              auto const& fzfab = flux[2]->array(););

    const Box& bxg2 = amrex::grow(bx,2);
    qtmp.resize(bxg2, NPRIM);
    auto const& q = qtmp.array();

    if (do_visc == 1)
    {
       diff_coeff.resize(bxg2, NCOEF);
    }

    amrex::ParallelFor(bxg2,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        cns_ctoprim(i, j, k, sfab, q, *lparm);
    });

    if (do_visc == 1)
    {
       auto const& coefs = diff_coeff.array();
       if(use_const_visc == 1 ) {
          amrex::ParallelFor(bxg2,
          [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
          {
              cns_constcoef(i, j, k, coefs, *lparm);
          });
       } else {
        //   amrex::ParallelFor(bxg2,
        //   [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        //   {
        //       cns_diffcoef(i, j, k, q, coefs, *lparm);
        //   });
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
    }

    const Box& bxg1 = amrex::grow(bx,1);
    slopetmp.resize(bxg1,UEINT);
    auto const& slope = slopetmp.array();

    auto l_plm_iorder = plm_iorder;
    auto l_plm_theta = plm_theta;

    // x-direction
    int cdir = 0;
    const Box& xslpbx = amrex::grow(bx, cdir, 1);
    amrex::ParallelFor(xslpbx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        cns_slope_x(i, j, k, slope, q, l_plm_iorder, l_plm_theta);
    });
    const Box& xflxbx = amrex::surroundingNodes(bx,cdir);
    amrex::ParallelFor(xflxbx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        cns_riemann_x(i, j, k, fxfab, slope, q, *lparm);
        for (int n = UEINT; n < NVAR; ++n) fxfab(i,j,k,n) = Real(0.0);
    });

    if (do_visc == 1)
    {
       auto const& coefs = diff_coeff.array();
       amrex::ParallelFor(xflxbx,
       [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
       {
           cns_diff_x(i, j, k, q, coefs, dxinv, fxfab);
       });
    }

    // y-direction
    cdir = 1;
    const Box& yslpbx = amrex::grow(bx, cdir, 1);
    amrex::ParallelFor(yslpbx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        cns_slope_y(i, j, k, slope, q, l_plm_iorder, l_plm_theta);
    });
    const Box& yflxbx = amrex::surroundingNodes(bx,cdir);
    amrex::ParallelFor(yflxbx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        cns_riemann_y(i, j, k, fyfab, slope, q, *lparm);
        for (int n = UEINT; n < NVAR; ++n) fyfab(i,j,k,n) = Real(0.0);
    });

    if(do_visc == 1)
    {
       auto const& coefs = diff_coeff.array();
       amrex::ParallelFor(yflxbx,
       [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
       {
           cns_diff_y(i, j, k, q, coefs, dxinv, fyfab);
       });
    }

#if (AMREX_SPACEDIM == 3)
    // z-direction
    cdir = 2;
    const Box& zslpbx = amrex::grow(bx, cdir, 1);
    amrex::ParallelFor(zslpbx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        cns_slope_z(i, j, k, slope, q, l_plm_iorder, l_plm_theta);
    });
    const Box& zflxbx = amrex::surroundingNodes(bx,cdir);
    amrex::ParallelFor(zflxbx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        cns_riemann_z(i, j, k, fzfab, slope, q, *lparm);
        for (int n = UEINT; n < NVAR; ++n) fzfab(i,j,k,n) = Real(0.0);
    });

    if(do_visc == 1)
    {
       auto const& coefs = diff_coeff.array();
       amrex::ParallelFor(zflxbx,
       [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
       {
           cns_diff_z(i, j, k, q, coefs, dxinv, fzfab);
       });
    }
#endif

    amrex::ParallelFor(bx, NVAR,
    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        cns_flux_to_dudt(i, j, k, n, dsdtfab, AMREX_D_DECL(fxfab,fyfab,fzfab), dxinv);
    });

//     if (gravity != Real(0.0)) {
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
