#include <AMReX_FluxRegister.H>
#include <CNS.h>
#include <prob.h>

#ifdef AMREX_USE_GPIBM
#include <ibm_solver.h>
#endif
using namespace amrex;

Real CNS::advance(Real time, Real dt, int /*iteration*/, int /*ncycle*/) {
  BL_PROFILE("CNS::advance()");

  state[0].allocOldData();
  state[0].swapTimeLevels(dt);
  
  MultiFab& S1 = get_old_data(State_Type);
  MultiFab& S2 = get_new_data(State_Type);

  int ncons = d_prob_closures->NCONS;
  int nghost= d_prob_closures->NGHOST;
  MultiFab Stemp(grids,dmap,ncons,nghost,MFInfo(),Factory());

  FluxRegister* fr_as_crse = nullptr;
  if (do_reflux && level < parent->finestLevel()) {
    CNS& fine_level = getLevel(level + 1);
    fr_as_crse = fine_level.flux_reg.get();
  }

  FluxRegister* fr_as_fine = nullptr;
  if (do_reflux && level > 0) {
    fr_as_fine = flux_reg.get();
  }

#ifdef AMREX_USE_GPIBM
  // Moving-geometry update.
  //
  // Geometry position is updated ONLY at the coarsest level (level 0).
  // The transform is global — it applies to all levels simultaneously.
  // Fine-level sub-cycles only rebuild markers and GPs at their own level
  // (the geometry position was already set by the coarse-level advance).
  //
  // This prevents fine-level sub-steps from advancing the geometry to
  // inconsistent times and avoids redundant transform updates.
  if (CNS::ib_move) {
    // Snapshot pre-move markers at THIS level
    auto& mfab_pre = *IBM::ib.bmf_a[level];
    FabArray<BaseFab<uint8_t>> old_markers(
        mfab_pre.boxArray(), mfab_pre.DistributionMap(),
        1, mfab_pre.nGrow(), MFInfo().SetArena(The_Managed_Arena()));
    for (MFIter mfi(mfab_pre, false); mfi.isValid(); ++mfi) {
      const Box& bx = mfi.fabbox();
      auto const& dst = old_markers.array(mfi);
      auto const& src = mfab_pre.const_array(mfi);
      ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        dst(i,j,k,0) = src(i,j,k,0);
      });
    }

    // Only the coarsest level updates the geometry transform.
    // Fine levels reuse the transform already set by the coarse advance.
    if (level == 0) {
#ifdef CNS_USE_FSI
#  ifdef CNS_FSI_DEFORMABLE
      PROB::update_geometry(time + dt, IBM::ib.geom_a, IBM::ib.ngeom);
      IBM::ib.rebuildGeometryData();
#  else
      PROB::update_rigid_transforms(time + dt, IBM::ib.transform_a, IBM::ib.ngeom);
      for (int i = 0; i < IBM::ib.ngeom; ++i) {
          IBM::ib.updateRigidTransform(i, IBM::ib.transform_a[i]);
      }
#  endif
#endif
    }

    // Rebuild markers and GPs at THIS level only (geometry is already at t^{n+1})
    rebuildIBM();

    // Fill cells newly exposed by the geometry motion at this level
    IBM::ib.fixExposedCells(old_markers, S1, level);
  }
#endif

  if (fr_as_crse) {
    fr_as_crse->setVal(Real(0.0));
  }

  if (order_rk == -2) {
    // Original time integration ///////////////////////////////////////////////
    // RK2 stage 1
    FillPatch(*this, Stemp, nghost, time, State_Type, 0, ncons);
    compute_rhs(Stemp, Real(0.5) * dt, fr_as_crse, fr_as_fine);
    // U^* = U^n + dt*dUdt^n
    MultiFab::LinComb(S2, Real(1.0), S1, 0, dt, Stemp, 0, 0, ncons, 0);
    // RK2 stage 2
    // After fillpatch Sborder = U^n+dt*dUdt^n
    state[0].setNewTimeLevel(time + dt);
    FillPatch(*this, Stemp, nghost, time + dt, State_Type, 0, ncons);
    compute_rhs(Stemp, Real(0.5) * dt, fr_as_crse, fr_as_fine);
    // S_new = 0.5*(Sborder+S_old) = U^n + 0.5*dt*dUdt^n
    MultiFab::LinComb(S2, Real(0.5), S1, 0, Real(0.5), S2, 0, 0, ncons, 0);
    // S_new += 0.5*dt*dSdt
    MultiFab::Saxpy(S2, Real(0.5) * dt, Stemp, 0, 0, ncons, 0);
    // We now have S_new = U^{n+1} = (U^n+0.5*dt*dUdt^n) + 0.5*dt*dUdt^*


    ////////////////////////////////////////////////////////////////////////////
  } else if (order_rk == 0) {  // returns rhs
    FillPatch(*this, Stemp, nghost, time, State_Type, 0, ncons);
    compute_rhs(Stemp, dt, fr_as_crse, fr_as_fine);
    MultiFab::Copy(S2, Stemp, 0, 0, ncons, 0);
  } else if (order_rk == 1) {
    FillPatch(*this, Stemp, nghost, time, State_Type, 0,
              ncons);  // filled at t_n to evalulate f(t_n,y_n).
    compute_rhs(Stemp, dt, fr_as_crse, fr_as_fine);
    MultiFab::LinComb(S2, Real(1.0), S1, 0, dt, Stemp, 0, 0, ncons, 0);
  } else if (order_rk == 2) {
    // Low-storage SSP-RK(m,2): m stages, C=m-1, C_eff=1-1/m.
    // Ref: Gottlieb et al., "Strong Stability Preserving Runge-Kutta and
    // Multistep Time Discretizations", §4.2.
    int m = stages_rk;
    MultiFab::Copy(S2, S1, 0, 0, ncons, 0);
    state[0].setOldTimeLevel(time);
    state[0].setNewTimeLevel(time);
    // First m-1 forward-Euler increments
    for (int i = 1; i <= m - 1; i++) {
      FillPatch(*this, Stemp, nghost, time + dt * Real(i - 1) / (m - 1),
                State_Type, 0, ncons);
      compute_rhs(Stemp, dt / Real(m - 1), fr_as_crse, fr_as_fine);
      MultiFab::Saxpy(S2, dt / Real(m - 1), Stemp, 0, 0, ncons, 0);
      state[State_Type].setNewTimeLevel(
          time + dt * Real(i) /
                     (m - 1));  // important to do this for correct fillpatch
                                // interpolations for the proceeding stages
    }
    // final stage
    FillPatch(*this, Stemp, nghost, time + dt, State_Type, 0, ncons);
    compute_rhs(Stemp, dt / Real(m - 1), fr_as_crse, fr_as_fine);
    MultiFab::LinComb(S2, Real(m - 1), S2, 0, dt, Stemp, 0, 0, ncons, 0);
    MultiFab::LinComb(S2, Real(1.0) / m, S1, 0, Real(1.0) / m, S2, 0, 0, ncons,
                      0);

    state[State_Type].setNewTimeLevel(time + dt);
  }

  // default
  else if (order_rk == 3) {
    if (stages_rk == 3) {
      // SSP-RK(3,3): http://ketch.github.io/numipedia/methods/SSPRK33.html
      state[0].setOldTimeLevel(time);
      FillPatch(*this, Stemp, nghost, time, State_Type, 0,
                ncons);  // filled at t_n to evalulate f(t_n,y_n).
      compute_rhs(Stemp, dt, fr_as_crse, fr_as_fine);
      MultiFab::LinComb(S2, Real(1.0), S1, 0, dt, Stemp, 0, 0, ncons, 0);

      state[0].setNewTimeLevel(
          time + dt);  // same time as upcoming FillPatch ensures we copy S2 to
                       // Sborder, without time interpolation
      FillPatch(*this, Stemp, nghost, time + dt, State_Type, 0, ncons);
      compute_rhs(Stemp, dt / 4, fr_as_crse, fr_as_fine);
      MultiFab::Xpay(Stemp, dt, S2, 0, 0, ncons, 0);
      MultiFab::LinComb(S2, Real(3.0) / 4, S1, 0, Real(1.0) / 4, Stemp, 0, 0,
                        ncons, 0);

      state[0].setNewTimeLevel(
          time + dt / 2);  // same time as upcoming FillPatch ensures we copy S2
                           // to Sborder, without time interpolation
      FillPatch(*this, Stemp, nghost, time + dt / 2, State_Type, 0, ncons);
      compute_rhs(Stemp, dt * Real(2.0) / 3, fr_as_crse, fr_as_fine);
      MultiFab::Xpay(Stemp, dt, S2, 0, 0, ncons, 0);
      MultiFab::LinComb(S2, Real(1.0) / 3, S1, 0, Real(2.0) / 3, Stemp, 0, 0,
                        ncons, 0);

      state[State_Type].setNewTimeLevel(
          time + dt);  // important to do this for correct fillpatch
                       // interpolations for the proceeding stages
    }

    else if (stages_rk == 4) {
      // SSP-RK(4,3): http://ketch.github.io/numipedia/methods/SSPRK43.html
      // Ref: Gottlieb et al., §4.2, p. 85.

      state[0].setOldTimeLevel(time);
      FillPatch(*this, Stemp, nghost, time, State_Type, 0, ncons);
      compute_rhs(Stemp, dt / 2, fr_as_crse, fr_as_fine);
      MultiFab::LinComb(S2, Real(1.0), S1, 0, dt / 2, Stemp, 0, 0, ncons, 0);

      state[0].setNewTimeLevel(
          time + dt / 2);  // same time as upcoming FillPatch ensures we copy S2
                           // to Sborder, without time interpolation
      FillPatch(*this, Stemp, nghost, time + dt / 2, State_Type, 0, ncons);
      compute_rhs(Stemp, dt / 2, fr_as_crse, fr_as_fine);
      MultiFab::Saxpy(S2, dt / 2, Stemp, 0, 0, ncons, 0);

      state[0].setNewTimeLevel(
          time + dt);  // same time as upcoming FillPatch ensures we copy S2 to
                       // Sborder, without time interpolation
      FillPatch(*this, Stemp, nghost, time + dt, State_Type, 0, ncons);
      compute_rhs(Stemp, dt / 6, fr_as_crse, fr_as_fine);
      MultiFab::LinComb(S2, Real(2.0) / 3, S1, 0, Real(1.0) / 3, S2, 0, 0,
                        ncons, 0);
      MultiFab::Saxpy(S2, dt / 6, Stemp, 0, 0, ncons, 0);

      state[0].setNewTimeLevel(
          time + dt / 2);  // same time as upcoming FillPatch ensures we copy S2
                           // to Sborder, without time interpolation
      FillPatch(*this, Stemp, nghost, time + dt / 2, State_Type, 0, ncons);
      compute_rhs(Stemp, dt / 2, fr_as_crse, fr_as_fine);
      MultiFab::Saxpy(S2, dt / 2, Stemp, 0, 0, ncons, 0);

      state[State_Type].setNewTimeLevel(
          time + dt);  // important to do this for correct fillpatch
                       // interpolations for the proceeding stages
    }

    else {
      // General SSP-RK(n^2, 3), n>2: C=2, C_eff=1-1/n (not yet implemented).
      // Ref: Gottlieb et al., §4.2, p. 85.
      amrex::Abort("SSPRK(n^2,3) with n>2 is not yet implemented");
    }

  }

#if ENSURE_MASSFRACSUM_ONE  
  clip_species_state(S2);
#endif

#ifdef AMREX_USE_GPIBM
  // ==========================================================================
  // End-of-step IBM correction (runs on S2 = state at t^{n+1})
  //
  // After the RK stages, the conservative state in IBM cells needs cleanup:
  //   [1] Ghost points: overwrite with wall-BC-reconstructed primitives
  //                     (always — also needed for static geometry).
  //   [2] Interior solid cells (FSI only): flood-fill from fluid/ghost
  //                     neighbors, then zero momentum. This keeps solid
  //                     cells carrying bounded, physically plausible data
  //                     so that AMR FillPatch/avgDown during regrid does
  //                     not interpolate garbage into fresh fluid cells.
  //   [3] Safety net:   replace any cell with catastrophically broken
  //                     values (NaN, Inf, or orders-of-magnitude outliers)
  //                     with a neighbor-based fallback. This should rarely
  //                     fire in a well-resolved simulation.
  // ==========================================================================
  {
    const PROB::ProbClosures& cls_h = *CNS::h_prob_closures;
    const PROB::ProbClosures* cls_d = CNS::d_prob_closures;
    auto& ib_mf = *IBM::ib.bmf_a[level];

    // ------------------------------------------------------------------------
    // Reconstruct ghost-point primitives from the current flow state + wall BC
    // ------------------------------------------------------------------------
    FillPatch(*this, Stemp, nghost, time + dt, State_Type, 0, ncons);

    MultiFab prims_mf(Stemp.boxArray(), Stemp.DistributionMap(),
                      cls_h.NPRIM, cls_h.NGHOST,
                      MFInfo().SetArena(The_Async_Arena()));
    for (MFIter mfi(Stemp, false); mfi.isValid(); ++mfi) {
      cls_h.cons2prims(mfi, Stemp.array(mfi), prims_mf.array(mfi));
    }

    // Sync wall-motion time so that compute_surfIB() sees the correct
    // instantaneous wall velocity for moving-wall BCs.
#ifdef CNS_USE_FSI
    PROB::Motion::sim_time = time + dt;
#endif

    IBM::ib.computeAllGPs(prims_mf, cls_d, level);
    Gpu::streamSynchronize();  // ensure GP primitives are fully written before Pass 1 reads them

    // ------------------------------------------------------------------------
    // Pass 1 (ALWAYS): Write GP-corrected primitives back to S2 as conservatives.
    //                   Only touches cells marked as ghost points (ibMarkers(,1)).
    // ------------------------------------------------------------------------
    for (MFIter mfi(S2, false); mfi.isValid(); ++mfi) {
      const Box& bx = mfi.tilebox();
      Array4<Real> const& state = S2.array(mfi);
      Array4<Real> const& prims = prims_mf.array(mfi);
      const auto& ibMarkers = ib_mf.array(mfi);

      ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        if (ibMarkers(i, j, k, 1)) {
          IntVect iv(AMREX_D_DECL(i, j, k));
          Real cons[PROB::ProbClosures::NCONS];
          cls_h.prims2cons(iv, prims, cons);
          for (int n = 0; n < PROB::ProbClosures::NCONS; n++) {
            state(i, j, k, n) = cons[n];
          }
        }
      });
    }

    // ------------------------------------------------------------------------
    // Pass 2 (ALWAYS when AMR is active): Flood-fill interior solid cells
    //                                     + zero their momentum.
    //
    // Why this is needed even for static geometry with AMR:
    //   Solid cells don't evolve in the RK step (RHS is zeroed), but WENO
    //   stencils near the surface still reach across the body and can
    //   contaminate interior solid cells with numerical garbage. These
    //   stale values then get read by avgDown/FillPatch during regrid and
    //   spread into fresh fluid cells on newly refined patches — which is
    //   exactly how AMR fails with T=0 outside the solid at step 1.
    //
    // For FSI (moving bodies), this is additionally critical because cells
    // that were solid can become fluid as the body moves.
    //
    // Zeroing momentum prevents spurious velocity amplification when
    // averaging fluid neighbors from different acoustic phases.
    // ------------------------------------------------------------------------
    if (parent->maxLevel() > 0 || CNS::ib_move) {
      for (MFIter mfi(S2, false); mfi.isValid(); ++mfi) {
        const Box& bx  = mfi.tilebox();
        const Box& bxg = mfi.growntilebox(d_prob_closures->NGHOST);
        auto const& state     = S2.array(mfi);
        auto const& ibMarkers = ib_mf.array(mfi);
        const int nc = ncons;

        // Tag array: 0 = already valid (fluid or GP-corrected ghost point)
        //            1 = interior solid, needs fixing
        //            2 = solid, already fixed in a previous iteration
        BaseFab<int> tagfab(bxg, 1, The_Managed_Arena());
        auto const& tag = tagfab.array();

        ParallelFor(bxg, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
          tag(i,j,k) = (ibMarkers(i,j,k,0) == 0 || ibMarkers(i,j,k,1) != 0)
                       ? 0 : 1;
        });

        // Iterative flood-fill: each iteration propagates valid data one
        // cell deeper into the solid region. Typically converges in 3-5
        // iterations; cap at 32 for safety.
        constexpr int MAX_FLOOD_ITER = 32;
        for (int iter = 0; iter < MAX_FLOOD_ITER; ++iter) {
          Gpu::DeviceScalar<int> d_nfixed(0);
          int* p_nfixed = d_nfixed.dataPtr();

          ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            if (tag(i,j,k) != 1) return;  // skip already-valid or already-fixed

            Real sum[PROB::ProbClosures::NCONS] = {};
            int count = 0;
            for (int dj = -1; dj <= 1; ++dj) {
              for (int di = -1; di <= 1; ++di) {
#if (AMREX_SPACEDIM == 3)
                for (int dk = -1; dk <= 1; ++dk) {
#else
                { int dk = 0;
#endif
                  if (di == 0 && dj == 0 && dk == 0) continue;
                  const int ii = i+di, jj = j+dj, kk = k+dk;
                  if (!bxg.contains(IntVect(AMREX_D_DECL(ii,jj,kk)))) continue;
                  if (tag(ii,jj,kk) == 0 || tag(ii,jj,kk) == 2) {
                    for (int n = 0; n < nc; ++n)
                      sum[n] += state(ii,jj,kk,n);
                    count++;
                  }
                }
              }
            }
            if (count > 0) {
              const Real inv = Real(1.0) / count;
              for (int n = 0; n < nc; ++n)
                state(i,j,k,n) = sum[n] * inv;
              tag(i,j,k) = 2;
              Gpu::Atomic::Add(p_nfixed, 1);
            }
          });

          Gpu::streamSynchronize();
          if (d_nfixed.dataValue() == 0) break;
        }

        // Zero momentum in interior solid cells.
        // Rationale: averaging fluid neighbors from different acoustic
        // phases can amplify velocity (observed solid |u| > fluid |u|
        // by 50-100%). During regrid, these spurious momenta feed into
        // FillPatch/avgDown and contaminate fresh fluid cells. Zeroing
        // is conservative and prevents this amplification loop.
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
          if (ibMarkers(i,j,k,0) == 0) return;  // fluid
          if (ibMarkers(i,j,k,1) != 0) return;  // ghost point (handled in Pass 1)

          using PC = PROB::ProbClosures;
          const Real rho = state(i,j,k, PC::URHO);
          if (rho <= Real(0)) return;

          // Subtract kinetic energy from total energy (keep internal energy)
          const Real mx = state(i,j,k, PC::UMX);
          const Real my = state(i,j,k, PC::UMY);
#if (AMREX_SPACEDIM == 3)
          const Real mz = state(i,j,k, PC::UMZ);
          const Real ke = Real(0.5) * (mx*mx + my*my + mz*mz) / rho;
#else
          const Real ke = Real(0.5) * (mx*mx + my*my) / rho;
#endif
          state(i,j,k, PC::UMX) = Real(0.0);
          state(i,j,k, PC::UMY) = Real(0.0);
#if (AMREX_SPACEDIM == 3)
          state(i,j,k, PC::UMZ) = Real(0.0);
#endif
          state(i,j,k, PC::UET) -= ke;
        });
      }
    } // end FSI flood-fill

    // ------------------------------------------------------------------------
    // Pass 3 (SAFETY NET, FSI only): Catch broken cells (NaN/Inf, non-positive,
    // or > 100x deviation from neighbor median) and replace with neighbor
    // median values + zero momentum.
    //
    // WHY FSI-ONLY: For static geometry, Pass 1 (GP correction) + Pass 2
    // (interior solid flood-fill) is sufficient. Unconditionally running
    // Pass 3 with its previous hardcoded magnitude bounds (rho in [0.01, 100],
    // rho_ref=1.177 SI air) silently replaces valid low-density cells in
    // non-dimensional cases (e.g. cylinder_Re40 with rho~0.003), causing
    // 300x density inflation and catastrophic eigenvalue blow-up. A
    // scale-adaptive magnitude check (deviation from neighbor median) is
    // also too aggressive in shock-dominated cases — early shock overshoots
    // trigger the check and get wrongly replaced, destroying the shock.
    //
    // For FSI, Pass 3 is still needed as a last-resort when cells transition
    // from solid to fluid with stale extrapolated values.
    if (CNS::ib_move) {
      constexpr Real MAG_TOL = Real(100.0);

      for (MFIter mfi(S2, false); mfi.isValid(); ++mfi) {
        const Box& bx  = mfi.tilebox();
        const Box& bxg = mfi.growntilebox(1);
        Array4<Real> const& state = S2.array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
          using PC = PROB::ProbClosures;
          const Real rho = state(i,j,k, PC::URHO);
          const Real E   = state(i,j,k, PC::UET);

          // Quick check: NaN/Inf/non-positive — always bad
          const bool hard_bad = !std::isfinite(rho) || rho <= Real(0.0)
                             || !std::isfinite(E)   || E   <= Real(0.0);

          // Collect valid neighbor densities and internal energies
          Real rho_nbrs[9] = {};
          Real eint_nbrs[9] = {};
          int n_nbrs = 0;
          for (int dj = -1; dj <= 1; ++dj) {
            for (int di = -1; di <= 1; ++di) {
#if (AMREX_SPACEDIM == 3)
              for (int dk = -1; dk <= 1; ++dk) {
#else
              { int dk = 0;
#endif
                if (di == 0 && dj == 0 && dk == 0) continue;
                if (n_nbrs >= 9) continue;
                const int ii = i+di, jj = j+dj, kk = k+dk;
                if (!bxg.contains(IntVect(AMREX_D_DECL(ii,jj,kk)))) continue;
                const Real rho_n = state(ii,jj,kk, PC::URHO);
                const Real E_n   = state(ii,jj,kk, PC::UET);
                if (!std::isfinite(rho_n) || rho_n <= Real(0.0)) continue;
                if (!std::isfinite(E_n)   || E_n   <= Real(0.0)) continue;
                const Real mx_n = state(ii,jj,kk, PC::UMX);
                const Real my_n = state(ii,jj,kk, PC::UMY);
#if (AMREX_SPACEDIM == 3)
                const Real mz_n = state(ii,jj,kk, PC::UMZ);
                const Real ke_n = Real(0.5) * (mx_n*mx_n + my_n*my_n + mz_n*mz_n) / rho_n;
#else
                const Real ke_n = Real(0.5) * (mx_n*mx_n + my_n*my_n) / rho_n;
#endif
                const Real eint_n = (E_n - ke_n) / rho_n;
                if (!std::isfinite(eint_n) || eint_n <= Real(0.0)) continue;
                rho_nbrs[n_nbrs]  = rho_n;
                eint_nbrs[n_nbrs] = eint_n;
                n_nbrs++;
              }
            }
          }

          if (n_nbrs == 0) return;  // no valid neighbors, leave alone

          // Simple median via partial sort
          for (int a = 0; a < n_nbrs - 1; ++a) {
            for (int b = a + 1; b < n_nbrs; ++b) {
              if (rho_nbrs[b] < rho_nbrs[a]) {
                Real t = rho_nbrs[a]; rho_nbrs[a] = rho_nbrs[b]; rho_nbrs[b] = t;
                t = eint_nbrs[a]; eint_nbrs[a] = eint_nbrs[b]; eint_nbrs[b] = t;
              }
            }
          }
          const Real rho_med  = rho_nbrs[n_nbrs / 2];
          const Real eint_med = eint_nbrs[n_nbrs / 2];

          // Scale-adaptive magnitude check
          bool mag_bad = false;
          if (!hard_bad) {
            const Real ratio = (rho > rho_med) ? rho / rho_med : rho_med / rho;
            if (ratio > MAG_TOL) mag_bad = true;
          }

          if (!hard_bad && !mag_bad) return;

          // Replace with neighbor-median values, zero momentum
          state(i,j,k, PC::URHO) = rho_med;
          state(i,j,k, PC::UMX)  = Real(0);
          state(i,j,k, PC::UMY)  = Real(0);
#if (AMREX_SPACEDIM == 3)
          state(i,j,k, PC::UMZ)  = Real(0);
#endif
          state(i,j,k, PC::UET)  = rho_med * eint_med;
        });
      }
    }
  }
#endif  // AMREX_USE_GPIBM

  return dt;
}
