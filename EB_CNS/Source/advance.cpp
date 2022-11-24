#if (AMREX_SPACEDIM > 1) //1D cannot have EB
#include <AMReX_EBFArrayBox.H>
#include <AMReX_MultiCutFab.H>
#endif

#include "CNS.H"

using namespace amrex;

/*
 * NOTE FOR FUTURE DEVELOPER(S): 
 *
 * When implementing stochastic fields, I use this file to contain all modifications.
 * The function advance deals with the entire state and reaction fab, while the 
 * compute_dSdt breaks the fab down into fields and feeds each field to the 
 * compute_dSdt_box(_eb) funstion. As a rule, only this file (and react.cpp) should  
 * contain stochastic field information, e.g. LEN_STATE, NUM_FIELD, and other files 
 * should see inputs of indivial field, e.g. NVAR, same as AMReX or PeleC.
 */

Real
CNS::advance(Real time, Real dt, int /*iteration*/, int /*ncycle*/)
{
    BL_PROFILE("CNS::advance()");

    // Prepare data fabs
    for (int i = 0; i < num_state_data_types; ++i) {
        if (!(i == Reactions_Type && do_react)) { // do not swap I_R (we only use new I_R data)
            state[i].allocOldData();
            state[i].swapTimeLevels(dt);
        }
    }

    MultiFab& S_new = get_new_data(State_Type);
    MultiFab& S_old = get_old_data(State_Type);
    MultiFab Sborder(grids, dmap, LEN_STATE, NUM_GROW, MFInfo(), Factory());
    MultiFab dSdt_old(grids, dmap, LEN_STATE, 0, MFInfo(), Factory());
    MultiFab dSdt_new(grids, dmap, LEN_STATE, 0, MFInfo(), Factory());
    dSdt_old.setVal(0.0);
    dSdt_new.setVal(0.0);

    MultiFab& I_R = get_new_data(Reactions_Type);

    get_new_data(Cost_Type).setVal(0.0);

    // Prepare flux register
#if (AMREX_SPACEDIM > 1) //1D cannot have EB
    EBFluxRegister* fr_as_crse = nullptr;
    EBFluxRegister* fr_as_fine = nullptr;
#else
    YAFluxRegister* fr_as_crse = nullptr;
    YAFluxRegister* fr_as_fine = nullptr;
#endif
    if (do_reflux && level < parent->finestLevel()) {
        fr_as_crse = &getLevel(level + 1).flux_reg;
    }
    if (do_reflux && level > 0) {
        fr_as_fine = &flux_reg;
    }

    if (fr_as_crse) {
        fr_as_crse->reset();
    }

    // Start time-stepping
    // RK2 stage 1: U^* = U^n + dt*dUdt^n + dt*I_R^n
    if (verbose > 0) amrex::Print() << " ... RK Step 1: Computing dSdt^{n}" << std::endl;

    FillPatch(*this, Sborder, NUM_GROW, time, State_Type, 0, LEN_STATE);
    // if (verbose >= 2) amrex::Print() << "                 ... RK Step 1.1" << std::endl;
    compute_dSdt(Sborder, dSdt_old, 0.5*dt, fr_as_crse, fr_as_fine, true);
    // if (verbose >= 2) amrex::Print() << "                 ... RK Step 1.2" << std::endl;
    MultiFab::LinComb(S_new, 1.0, Sborder, 0, dt, dSdt_old, 0, 0, LEN_STATE, 0);
    // if (verbose >= 2) amrex::Print() << "                 ... RK Step 1.3" << std::endl;
    if (do_react) {
        MultiFab::Saxpy(S_new, dt, I_R, 0, UFS, NUM_SPECIES, 0);
        MultiFab::Saxpy(S_new, dt, I_R, NUM_SPECIES, UEDEN, 1, 0);
    }
    computeTemp(S_new, 0); // Update EINT and TEMP
    
    // RK2 stage 2: U^{n+1} = U^n + 0.5*dt*(dUdt^n + dUdt^{n+1}) + dt*I_R^{n+1}
    if (verbose > 0) amrex::Print() << " ... RK Step 2: Computing dSdt^{n+1}" << std::endl;

    FillPatch(*this, Sborder, NUM_GROW, time+dt, State_Type, 0, LEN_STATE);
    compute_dSdt(Sborder, dSdt_new, 0.5*dt, fr_as_crse, fr_as_fine, (rk_reaction_iter < 1));
    MultiFab::LinComb(S_new, 1.0, S_old, 0, 0.5*dt, dSdt_old, 0, 0, LEN_STATE, 0);
    MultiFab::Saxpy(S_new, 0.5*dt, dSdt_new, 0, 0, LEN_STATE, 0); // U^** = U^n + 0.5*dt*(dUdt^n + dUdt^{n+1})

    if (do_react) {
        react_source(time, dt, false); // Compute I_R^{n+1}(U^**)

        // Finally, U^{n+1} = U^** + dt*I_R^{n+1}
        MultiFab::Saxpy(S_new, dt, I_R, 0, UFS, NUM_SPECIES, 0);
        MultiFab::Saxpy(S_new, dt, I_R, NUM_SPECIES, UEDEN, 1, 0);
    }
    computeTemp(S_new, 0); // Update EINT and TEMP
    
    // Iterate to couple chemistry
    if (do_react && (rk_reaction_iter > 1)) {
        for (int iter = 1; iter < rk_reaction_iter; ++iter) {
            if (verbose > 0) {
                amrex::Print() << " ... Re-computing dSdt^{n+1}: iter " 
                               << iter << " of " << rk_reaction_iter << ")" << std::endl;
            }

            FillPatch(*this, Sborder, NUM_GROW, time+dt, State_Type, 0, LEN_STATE);
            compute_dSdt(Sborder, dSdt_new, 0.5*dt, fr_as_crse, fr_as_fine, (iter == rk_reaction_iter));
            
            MultiFab::LinComb(S_new, 1.0, S_old, 0, 0.5*dt, dSdt_old, 0, 0, LEN_STATE, 0);
            MultiFab::Saxpy(S_new, 0.5*dt, dSdt_new, 0, 0, LEN_STATE, 0);
            
            react_source(time, dt, false);

            MultiFab::Saxpy(S_new, dt, I_R, 0, UFS, NUM_SPECIES, 0);
            MultiFab::Saxpy(S_new, dt, I_R, NUM_SPECIES, UEDEN, 1, 0);
            
            computeTemp(S_new, 0);
        }
    }

    return dt;
}

void
CNS::compute_dSdt (const MultiFab& S, MultiFab& dSdt, Real dt,
#if (AMREX_SPACEDIM > 1) //1D cannot have EB
                   EBFluxRegister* fr_as_crse, EBFluxRegister* fr_as_fine,
#else
                   YAFluxRegister* fr_as_crse, YAFluxRegister* fr_as_fine,
#endif
                   bool write_to_flux_register)
{
    BL_PROFILE("CNS::compute_dSdt()");

    const Real* dx = geom.CellSize();
    const int ncomp = dSdt.nComp();

    int as_crse = (fr_as_crse != nullptr);
    int as_fine = (fr_as_fine != nullptr);

    MultiFab& cost = get_new_data(Cost_Type);

#if (AMREX_SPACEDIM > 1) //1D cannot have EB
    auto const& fact = dynamic_cast<EBFArrayBoxFactory const&>(S.Factory());
    auto const& flags = fact.getMultiEBCellFlagFab();
#endif

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
{
    std::array<FArrayBox,AMREX_SPACEDIM> flux;
    FArrayBox dm_as_fine(Box::TheUnitBox(), ncomp);
    FArrayBox fab_drho_as_crse(Box::TheUnitBox(), ncomp);
    IArrayBox fab_rrflag_as_crse(Box::TheUnitBox());
    
    for (MFIter mfi(S, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        auto wt = amrex::second();

        const Box& bx = mfi.tilebox();
#if (AMREX_SPACEDIM > 1) //1D cannot have EB
        const auto& flag = flags[mfi];
        if (flag.getType(bx) == FabType::covered) {
            dSdt[mfi].setVal<RunOn::Device>(0.0, bx, 0, ncomp);
        } 
        else 
#endif
        {
            // flux is used to store centroid flux needed for reflux
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                flux[idim].resize(amrex::surroundingNodes(bx,idim), LEN_STATE);
                flux[idim].setVal<RunOn::Device>(0.);
            }

#if (AMREX_SPACEDIM > 1) //1D cannot have EB
            if (flag.getType(amrex::grow(bx,1)) == FabType::regular) {
#endif
                Array4<Real const>    s_arr =    S.array(mfi);
                Array4<Real      > dsdt_arr = dSdt.array(mfi);
                // AMREX_D_TERM(
                //     Array4<Real  > xflx_arr = flux[0].array(); ,
                //     Array4<Real  > yflx_arr = flux[1].array(); ,
                //     Array4<Real  > zflx_arr = flux[2].array();)

                compute_dSdt_box(bx, s_arr, dsdt_arr, {AMREX_D_DECL(&flux[0],&flux[1],&flux[2])});

                if (write_to_flux_register) {
                    if (fr_as_crse) {
                        fr_as_crse->CrseAdd(mfi, {AMREX_D_DECL(&flux[0], &flux[1], &flux[2])}, dx, dt, RunOn::Device);
                    }

                    if (fr_as_fine) {
                        fr_as_fine->FineAdd(mfi, {AMREX_D_DECL(&flux[0], &flux[1], &flux[2])}, dx, dt, RunOn::Device);
                    }
                }
#if (AMREX_SPACEDIM > 1) //1D cannot have EB
            } else {
                FArrayBox* p_drho_as_crse = (fr_as_crse) ?
                    fr_as_crse->getCrseData(mfi) : &fab_drho_as_crse;
                const IArrayBox* p_rrflag_as_crse = (fr_as_crse) ?
                    fr_as_crse->getCrseFlag(mfi) : &fab_rrflag_as_crse;

                if (fr_as_fine) {
                    dm_as_fine.resize(amrex::grow(bx,1),ncomp);
                }

                Array4<Real const> vf_arr = (*volfrac).array(mfi);
                Array4<Real const> bcent_arr = (*bndrycent).array(mfi);

                AMREX_D_TERM(Array4<Real const> const& apx = areafrac[0]->const_array(mfi);,
                             Array4<Real const> const& apy = areafrac[1]->const_array(mfi);,
                             Array4<Real const> const& apz = areafrac[2]->const_array(mfi));
                AMREX_D_TERM(Array4<Real const> const& fcx = facecent[0]->const_array(mfi);,
                             Array4<Real const> const& fcy = facecent[1]->const_array(mfi);,
                             Array4<Real const> const& fcz = facecent[2]->const_array(mfi));
                
                Array4<Real const> const&    s_arr =    S.array(mfi);
                Array4<Real      > const& dsdt_arr = dSdt.array(mfi);
                AMREX_D_TERM(Array4<Real> xflx_arr = flux[0].array(); ,
                             Array4<Real> yflx_arr = flux[1].array(); ,
                             Array4<Real> zflx_arr = flux[2].array();)
                
                compute_dSdt_box_eb(bx, s_arr, dsdt_arr,
                                    AMREX_D_DECL(xflx_arr, yflx_arr, zflx_arr),
                                    flags.const_array(mfi), vf_arr,
                                    AMREX_D_DECL(apx, apy, apz), AMREX_D_DECL(fcx, fcy, fcz), bcent_arr,
                                    as_crse, p_drho_as_crse->array(), p_rrflag_as_crse->const_array(),
                                    as_fine, dm_as_fine.array(), level_mask.const_array(mfi), dt);
                
                // // dSdt_Fields
                // for (int nf = 0; nf < NUM_FIELD; ++nf) {
                //     Array4<Real const> const&    s_arr =    S[mfi].array(NVAR + nf*NVAR);
                //     Array4<Real      > const& dsdt_arr = dSdt[mfi].array(NVAR + nf*NVAR);
                //     AMREX_D_TERM(
                //         Array4<Real  > xflx_arr = flux[0].array(NVAR + nf*NVAR); ,
                //         Array4<Real  > yflx_arr = flux[1].array(NVAR + nf*NVAR); ,
                //         Array4<Real  > zflx_arr = flux[2].array(NVAR + nf*NVAR);)
                //     compute_dSdt_box_eb(bx, s_arr, dsdt_arr,
                //                         AMREX_D_DECL(xflx_arr, yflx_arr, zflx_arr),
                //                         flags.const_array(mfi), vf_arr,
                //                         AMREX_D_DECL(apx,apy,apz),AMREX_D_DECL(fcx,fcy,fcz), bcent_arr,
                //                         as_crse, p_drho_as_crse->array(NVAR + nf*NVAR), p_rrflag_as_crse->const_array(),
                //                         as_fine, dm_as_fine.array(NVAR + nf*NVAR), level_mask.const_array(mfi), dt);
                // }

                if (write_to_flux_register) {
                    if (fr_as_crse) {
                        fr_as_crse->CrseAdd(mfi, {AMREX_D_DECL(&flux[0], &flux[1], &flux[2])},
                                            dx, dt, (*volfrac)[mfi],
                                            {AMREX_D_DECL(&((*areafrac[0])[mfi]),
                                                          &((*areafrac[1])[mfi]),
                                                          &((*areafrac[2])[mfi]))},
                                            RunOn::Device);
                    }

                    if (fr_as_fine) {
                        fr_as_fine->FineAdd(mfi, {AMREX_D_DECL(&flux[0], &flux[1], &flux[2])},
                                            dx, dt, (*volfrac)[mfi],
                                            {AMREX_D_DECL(&((*areafrac[0])[mfi]),
                                                          &((*areafrac[1])[mfi]),
                                                          &((*areafrac[2])[mfi]))},
                                            dm_as_fine,
                                            RunOn::Device);
                    }
                }
            }
#endif
        }
        
        Gpu::streamSynchronize();

        wt = (amrex::second() - wt) / bx.d_numPts();
        cost[mfi].plus<RunOn::Device>(wt, bx);
    }
} // omp parallel
}
