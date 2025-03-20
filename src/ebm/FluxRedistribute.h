#ifndef FluxRedistribute_H_
#define FluxRedistribute_H_

#include <AMReX_EBFluxRegister.H>
#include <AMReX_YAFluxRegister.H>
#include <AMReX_EBFArrayBox.H>
#include <AMReX_MultiCutFab.H>
#include <AMReX_YAFluxRegister_K.H>


inline void cerisse_flux_redistribute (
    const Box& bx,
    Array4<Real            > const& dqdt,
    Array4<Real       const> const& divc,
    Array4<Real       const> const& wt,
    Array4<Real       const> const& vfrac,
    Array4<EBCellFlag const> const& flag,
    const Geometry& geom,
    int ncomp, Real dt)
{
    //
    // Assume grid is uniform
    //
    const Real* dx = geom.CellSize();

    const Box dbox1 = geom.growPeriodicDomain(1);
    const Box dbox2 = geom.growPeriodicDomain(2);

    const Box& grown1_bx = amrex::grow(bx,1);
    const Box& grown2_bx = amrex::grow(bx,2);

    Real reredistribution_threshold = amrex_eb_get_reredistribution_threshold();

    int bx_ilo = bx.smallEnd()[0];
    int bx_ihi = bx.bigEnd()[0];
    int bx_jlo = bx.smallEnd()[1];
    int bx_jhi = bx.bigEnd()[1];
#if (AMREX_SPACEDIM == 3)
    int bx_klo = bx.smallEnd()[2];
    int bx_khi = bx.bigEnd()[2];
#endif

    //
    // Working arrays
    //
    FArrayBox  delm_fab(grown1_bx,ncomp);
    FArrayBox  optmp_fab(grown2_bx,ncomp);
    FArrayBox  mask_fab(grown2_bx);
    FArrayBox  fluid_fab(grown2_bx);
    

    Array4<Real> const& optmp = optmp_fab.array();
    Array4<Real> const& mask  = mask_fab.array();
    Array4<Real> const& fluid = fluid_fab.array();

    Array4<Real> const& delm  = delm_fab.array();

    //
    // Array "mask" is used to sever the link to ghost cells when the BCs
    // are not periodic
    // mask  is  1 when a cell can be used in computations, 0 otherwise
    // fluid is  1 when a cell is not covered, 0 otherwise

    //
    AMREX_FOR_3D(grown2_bx, i, j, k,
    {
      mask(i,j,k) = (dbox2.contains(IntVect(AMREX_D_DECL(i,j,k)))) ? 1.0 : 0.0;
      fluid(i,j,k) =  !flag(i,j,k).isCovered() ? 1.0: 0.0;
    });

    //
    // Init to zero  DRHS array
    //
    AMREX_FOR_4D(grown2_bx, ncomp, i, j, k, n,
    {
      optmp(i,j,k,n) = 0;
    });

    //const Real *prob_lo = geomdata.ProbLo();

    //
    // Step 2: compute delta M (mass gain or loss) on (lo-1,lo+1) grown1_bx
    //
      
    
    constexpr int nb = 1;
    constexpr Real epsdx = 0.01;

    amrex::ParallelFor(grown1_bx, ncomp,
    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
      // check if it is cut cell
      if (flag(i,j,k).isSingleValued())
      {
        Real vtot(0.);
        Real divnc(0.);

        // Real x = prob_lo[0] + (i + Real(0.5)) * dx[0];
        // Real y = prob_lo[1] + (j + Real(0.5)) * dx[1];

#if (AMREX_SPACEDIM == 2)
        int kk(0);
#else
        for (int kk = -nb; kk <= nb; kk++) {
#endif
        for (int jj = -nb; jj <= nb; jj++) {
        for (int ii = -nb; ii <= nb; ii++) {
          // loop over neighbours
          // check connected neighbours except itself  and the box contains node

        //  if ( (ii != 0 || jj != 0 || kk != 0) && flag(i,j,k).isConnected(ii,jj,kk) &&
        //                  dbox2.contains(IntVect(AMREX_D_DECL(i+ii,j+jj,k+kk))))

  //        if (flag(i,j,k).isConnected(ii,jj,kk) )            
            {
              Real drx  = ii*Real(0.5)* dx[0];
              Real dry  = jj*Real(0.5)* dx[1];
#if (AMREX_SPACEDIM == 2)
              Real drz  = 0;
#else
              Real drz  = kk*Real(0.5)* dx[2];              
#endif

              Real wcorr =  1.0/(epsdx*dx[0] + std::sqrt(drx*drx + dry*dry + drz*drz));
              Real unwted_frac = vfrac(i+ii,j+jj,k+kk) * mask(i+ii,j+jj,k+kk)*wcorr;
              vtot  += unwted_frac;
              divnc += unwted_frac*divc(i+ii,j+jj,k+kk,n);
            }

        AMREX_D_TERM(},},}) // end loop nb
        divnc /= vtot;
        // We need to multiply by mask to make sure optmp is zero for cells
        // outside the domain for non-cyclic BCs
        optmp(i,j,k,n) =  (1 - vfrac(i,j,k)) * (divnc - divc(i,j,k,n)) * mask(i,j,k);
        delm(i,j,k,n)  = -(    vfrac(i,j,k)) * optmp(i,j,k,n);
      } else {
        delm(i,j,k,n) = 0.;
      }
    });
    

    //
    // Step 3: redistribute excess/loss of mass
    //
    //

    amrex::ParallelFor(grown1_bx, ncomp,
    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
      bool valid_dst_cell;
       // check if it is cut cell
      if (flag(i,j,k).isSingleValued())
        {
          Real wtot = Real(0.);
          // calculate w  = sum_nb volfrac(nb)*wt(nb) 
#if (AMREX_SPACEDIM == 2)
          int kk(0);
#else
          for (int kk = -nb; kk <= nb; kk++) {
#endif
          for (int jj = -nb; jj <= nb; jj++) {
          for (int ii = -nb; ii <= nb; ii++) {
          // loop over neighbours  
          // check connected neighbours except itself  
          // if ( (ii != 0 || jj != 0 || kk != 0) && flag(i,j,k).isConnected(ii,jj,kk) )

          //  if (flag(i,j,k).isConnected(ii,jj,kk) ) 
              {
                wtot += vfrac(i+ii,j+jj,k+kk)*wt(i+ii,j+jj,k+kk)* mask(i+ii,j+jj,k+kk);
              }
          AMREX_D_TERM(},},})

           // calculate 1/w  = 1/sum_nb volfrac(nb)*wt(nb) 
#ifdef AMREX_USE_FLOAT
            wtot = Real(1.0)/(wtot + Real(1.e-30));
#else
            wtot = Real(1.0)/(wtot + Real(1.e-80));
#endif


#if (AMREX_SPACEDIM == 2)
            kk = 0;
#else
            for (int kk = -nb; kk <= nb; kk++) {
#endif
            for (int jj = -nb; jj <= nb; jj++) {
            for (int ii = -nb; ii <= nb; ii++) {
              // loop over neighbours  
              // check connected neighbours except itself  
              //  if ( (ii != 0 || jj != 0 || kk != 0) && flag(i,j,k).isConnected(ii,jj,kk) )
              // if (flag(i,j,k).isConnected(ii,jj,kk) ) 
                {
                  int iii = i + ii;
                  int jjj = j + jj;
                  int kkk = k + kk;

                  Real drho = delm(i,j,k,n)*wtot*wt(iii,jjj,kkk)* mask(iii,jjj,kkk) ;

                //    Gpu::Atomic::Add(&optmp(iii,jjj,kkk,n), drho); // SNM

                   optmp(iii,jjj,kkk,n) +=drho;

                } // isConnected
            AMREX_D_TERM(},},})
        } // isSingleValued
    });

    // RHS = RHS(uncorrected) + DELTA_RHS
    amrex::ParallelFor(bx, ncomp,
    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        // if (!flag(i,j,k).isCovered()) {
        //     dqdt(i,j,k,n) = divc(i,j,k,n) + optmp(i,j,k,n);
        // }

        dqdt(i,j,k,n) = divc(i,j,k,n) + optmp(i,j,k,n)*fluid(i,j,k);

    });


    Gpu::streamSynchronize(); // because of FArrayBoxes defined in this function
} // end cerisse_flux_redistribute

//
// Do small cell redistribution on one FAB with the Array4's already passed in
//

    // int as_crse = 0;
    // int as_fine = 0;
    // Real dummy_dt  = 0.0;
    // Array4<Real> dummy_drho_crse  = Array4<Real>();
    // Array4<Real> dummy_dm_as_fine = Array4<Real>();
    // Array4<int > dummy_levmsk     = Array4<int>();
    // Array4<int > dummy_flag_crse  = Array4<int>();
    // int dummy_not_covered = 0;
    // amrex_flux_redistribute (bx, div, divc, wt, vfrac, flag_arr,
    //                          as_crse, dummy_drho_crse , dummy_flag_crse,
    //                          as_fine, dummy_dm_as_fine, dummy_levmsk,
    //                          geom, use_wts_in_divnc, dummy_not_covered,
    //                          icomp, ncomp, dummy_dt);
// end namespace

#endif
