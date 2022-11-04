#include "CNS.H"

#include <AMReX_EBFArrayBox.H>
#include <AMReX_MultiCutFab.H>

using namespace amrex;

Real
CNS::advance(Real time, Real dt, int /*iteration*/, int /*ncycle*/)
{
    BL_PROFILE("CNS::advance()");

    for (int i = 0; i < num_state_data_types; ++i) {
        if (!(i == Reactions_Type && do_react)) { // do not swap I_R yet
            state[i].allocOldData();
            state[i].swapTimeLevels(dt);
        }
    }

    MultiFab& S_new = get_new_data(State_Type);
    MultiFab& S_old = get_old_data(State_Type);
    MultiFab dSdt(grids, dmap, LEN_STATE, 0, MFInfo(), Factory());
    MultiFab Sborder(grids, dmap, LEN_STATE, NUM_GROW, MFInfo(), Factory());
    dSdt.setVal(0.);

    MultiFab& I_R = get_new_data(Reactions_Type);

    get_new_data(Cost_Type).setVal(0.0);

    EBFluxRegister* fr_as_crse = nullptr;
    if (do_reflux && level < parent->finestLevel()) {
        fr_as_crse = &getLevel(level + 1).flux_reg;
    }

    EBFluxRegister* fr_as_fine = nullptr;
    if (do_reflux && level > 0) {
        fr_as_fine = &flux_reg;
    }

    if (fr_as_crse) {
        fr_as_crse->reset();
    }

    // RK2 stage 1
    // Sborder = S_old = U^n
    FillPatch(*this, Sborder, NUM_GROW, time, State_Type, 0, LEN_STATE);
    compute_dSdt(Sborder, dSdt, 0.5*dt, fr_as_crse, fr_as_fine);
                             // ^ should this be 0.5*dt or dt?                             
    // U^* = U^n + dt*dUdt^n + dt*I_R^{n-1/2}
    MultiFab::LinComb(S_new, 1.0, Sborder, 0, dt, dSdt, 0, 0, LEN_STATE, 0);
    if (do_react) {
        for (int nf = 0; nf < NUM_FIELD; ++nf) {// src start         des start   len          ng
            MultiFab::Saxpy(S_new, dt, I_R, nf * NREACT,      nf * UFS,   NUM_SPECIES, 0);
            MultiFab::Saxpy(S_new, dt, I_R, nf * NUM_SPECIES, nf * UEDEN, 1,           0);
        }        
    }
    computeTemp(S_new, 0);

    // RK2 stage 2
    // Sborder = U^*
    FillPatch(*this, Sborder, NUM_GROW, time+dt, State_Type, 0, LEN_STATE);
    compute_dSdt(Sborder, dSdt, 0.5*dt, fr_as_crse, fr_as_fine);
                             // ^ should this be 0.5*dt or dt?
    // U^** = U^n + 0.5*dt*(dUdt^n + dUdt^{n+1}) + 0.5*dt*I_R^{n-1/2} = 0.5*(Sborder+S_old) + 0.5*dt*dSdt
    MultiFab::LinComb(S_new, 0.5, Sborder, 0, 0.5, S_old, 0, 0, LEN_STATE, 0);
    MultiFab::Saxpy(S_new, 0.5*dt, dSdt, 0, 0, LEN_STATE, 0);

    if (do_react) {
        // U^** = U^n + 0.5*dt*(dUdt^n + dUdt^{n+1}) + dt*I_R^{n-1/2}
        for (int nf = 0; nf < NUM_FIELD; ++nf) {
            MultiFab::Saxpy(S_new, 0.5 * dt, I_R, nf * NREACT, nf * UFS, NUM_SPECIES, 0);
            MultiFab::Saxpy(S_new, 0.5 * dt, I_R, nf * NUM_SPECIES, nf * UEDEN, 1, 0);
        }

        // dSdt_AD = (1/dt)(U^** - U^n) - I_R^n = 0.5*(dUdt^n + dUdt^{n+1})
        MultiFab::LinComb(dSdt, 1.0/dt, S_new, 0, -1.0/dt, S_old, 0, 0, LEN_STATE, 0);
        for (int nf = 0; nf < NUM_FIELD; ++nf) {
            MultiFab::Subtract(dSdt, I_R, nf * NREACT, nf * UFS, NUM_SPECIES, 0);
            MultiFab::Subtract(dSdt, I_R, nf * NUM_SPECIES, nf * UEDEN, 1, 0);
        }

        // Compute I_R^{n+1/2}
        reaction_source(I_R, S_new, time, dt, false, &dSdt);

        // Finally, U^{n+1} = U^n + dt*(dSdt_AD + I_R^{n+1}) = U^** + dt*I_R^{n+1/2}
        MultiFab::LinComb(S_new, 1.0, S_old, 0, dt, dSdt, 0, 0, NVAR, 0);
        for (int nf = 0; nf < NUM_FIELD; ++nf) {
            MultiFab::Saxpy(S_new, dt, I_R, nf * NREACT, nf * UFS, NUM_SPECIES, 0);
            MultiFab::Saxpy(S_new, dt, I_R, nf * NUM_SPECIES, nf * UEDEN, 1, 0);
        }
    }
    computeTemp(S_new, 0);    

    return dt;
}

void
CNS::compute_dSdt(const MultiFab& S, MultiFab& dSdt, Real dt,
                  EBFluxRegister* fr_as_crse, EBFluxRegister* fr_as_fine)
{
    BL_PROFILE("CNS::compute_dSdt()");

    const Real* dx = geom.CellSize(); //dir? CellSizeArray()?
    const int ncomp = dSdt.nComp();

    int as_crse = (fr_as_crse != nullptr);
    int as_fine = (fr_as_fine != nullptr);

    MultiFab& cost = get_new_data(Cost_Type);

    auto const& fact = dynamic_cast<EBFArrayBoxFactory const&>(S.Factory());
    auto const& flags = fact.getMultiEBCellFlagFab();

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

            const auto& flag = flags[mfi];

            if (flag.getType(bx) == FabType::covered) {
                dSdt[mfi].setVal<RunOn::Device>(0.0, bx, 0, ncomp);
            } else {
                // flux is used to store centroid flux needed for reflux
                for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                    flux[idim].resize(amrex::surroundingNodes(bx,idim), LEN_STATE);
                    flux[idim].setVal<RunOn::Device>(0.);
                }

                if (flag.getType(amrex::grow(bx,1)) == FabType::regular)
                {
                    for (int nf = 0; nf < NUM_FIELD; ++nf) {
                        Array4<Real const>    s_arr =    S[mfi].array(nf * NVAR);
                        Array4<Real      > dsdt_arr = dSdt[mfi].array(nf * NVAR);
                        AMREX_D_TERM(
                            Array4<Real  > xflx_arr = flux[0].array(nf * NVAR); ,
                            Array4<Real  > yflx_arr = flux[1].array(nf * NVAR); ,
                            Array4<Real  > zflx_arr = flux[2].array(nf * NVAR);)
                    
                        compute_dSdt_box(bx, s_arr, dsdt_arr, AMREX_D_DECL(xflx_arr, yflx_arr, zflx_arr));
                    }

                    if (fr_as_crse) {
                        fr_as_crse->CrseAdd(mfi, {AMREX_D_DECL(&flux[0], &flux[1], &flux[2])}, dx, dt, RunOn::Device);
                    }

                    if (fr_as_fine) {
                        fr_as_fine->FineAdd(mfi, {AMREX_D_DECL(&flux[0], &flux[1], &flux[2])}, dx, dt, RunOn::Device);
                    }
                }
                else
                {
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
                    
                    for (int nf = 0; nf < NUM_FIELD; ++nf) {
                        Array4<Real const> const&    s_arr =    S[mfi].array(nf * NVAR);
                        Array4<Real      > const& dsdt_arr = dSdt[mfi].array(nf * NVAR);
                        AMREX_D_TERM(
                            Array4<Real  > xflx_arr = flux[0].array(nf * NVAR); ,
                            Array4<Real  > yflx_arr = flux[1].array(nf * NVAR); ,
                            Array4<Real  > zflx_arr = flux[2].array(nf * NVAR);)

                        compute_dSdt_box_eb(bx, s_arr, dsdt_arr,
                                            AMREX_D_DECL(xflx_arr, yflx_arr, zflx_arr),
                                            flags.const_array(mfi), vf_arr,
                                            AMREX_D_DECL(apx,apy,apz),AMREX_D_DECL(fcx,fcy,fcz), bcent_arr,
                                            as_crse, p_drho_as_crse->array(nf * NVAR), p_rrflag_as_crse->const_array(),
                                            as_fine, dm_as_fine.array(nf * NVAR), level_mask.const_array(mfi), dt);
                    }

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

            Gpu::streamSynchronize();

            wt = (amrex::second() - wt) / bx.d_numPts();
            cost[mfi].plus<RunOn::Device>(wt, bx);
        }
    }
}
