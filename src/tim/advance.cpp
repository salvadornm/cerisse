#include <AMReX_FluxRegister.H>
#include <CNS.h>
#include <CNSconstants.h>
#include <prob.h>
#include <iomanip>
// 
#include "mandebug.h"

#ifdef AMREX_USE_GPIBM
#include <ibm_solver.h>
#endif
using namespace amrex;

namespace {

using FaceData = CNS::NSCBCFaceData;

void define_face_data_like( FaceData& destination, FaceData const& source, int ncomp, int ngrow = 0)
{
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {

        if (!source[dir]) {
            destination[dir].reset();
            continue;
        }

        destination[dir] =
            std::make_unique<amrex::MultiFab>(
                source[dir]->boxArray(),
                source[dir]->DistributionMap(),
                ncomp,
                ngrow);

        destination[dir]->setVal(amrex::Real(0.0));
    }
}

void copy_face_data(
    FaceData& destination,
    FaceData const& source,
    int ncomp)
{
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {

        if (!destination[dir] || !source[dir]) {
            continue;
        }

        amrex::MultiFab::Copy(
            *destination[dir],
            *source[dir],
            0,
            0,
            ncomp,
            0);
    }
}

void saxpy_face_data(
    FaceData& destination,
    amrex::Real factor,
    FaceData const& source,
    int ncomp)
{
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {

        if (!destination[dir] || !source[dir]) {
            continue;
        }

        amrex::MultiFab::Saxpy(
            *destination[dir],
            factor,
            *source[dir],
            0,
            0,
            ncomp,
            0);
    }
}

void lincomb_face_data(
    FaceData& destination,
    amrex::Real factor_a,
    FaceData const& source_a,
    amrex::Real factor_b,
    FaceData const& source_b,
    int ncomp)
{
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {

        if (!destination[dir] ||
            !source_a[dir] ||
            !source_b[dir]) {
            continue;
        }

        amrex::MultiFab::LinComb(
            *destination[dir],
            factor_a,
            *source_a[dir],
            0,
            factor_b,
            *source_b[dir],
            0,
            0,
            ncomp,
            0);
    }
}

} // namespace


Real CNS::advance(Real time, Real dt, int /*iteration*/, int /*ncycle*/) {
  BL_PROFILE("CNS::advance()");

  state[0].allocOldData();
  state[0].swapTimeLevels(dt);
  
  MultiFab& S1 = get_old_data(State_Type);
  MultiFab& S2 = get_new_data(State_Type);

  int ncons = d_prob_closures->NCONS;
  int nghost= d_prob_closures->NGHOST;

  // temp
  // int nglin = 0; int ngfill = nghost;
  // if (CNS::use_nscbc) { nglin = nghost; ngfill = 0;}


  MultiFab Stemp(grids,dmap,ncons,nghost,MFInfo(),Factory());

  FluxReg* fr_as_crse = nullptr;
  if (do_reflux && level < parent->finestLevel()) {
    fr_as_crse = &getLevel(level + 1).flux_reg;
  }

  FluxReg* fr_as_fine = nullptr;
  if (do_reflux && level > 0) {
    fr_as_fine = &flux_reg;
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
    fr_as_crse->reset();
  }



  ////////////////////////////////////////////////////////////////////////////
  if (order_rk == -2) {
    
    // Temporary arrays
    NSCBCFaceData UBC_old; NSCBCFaceData UBC_rhs1; NSCBCFaceData UBC_rhs2;

    if (use_nscbc) {
      define_face_data_like( UBC_old, nscbc_ubc, ncons);
      define_face_data_like( UBC_rhs1,nscbc_ubc, ncons);
      define_face_data_like( UBC_rhs2,nscbc_ubc, ncons);
      copy_face_data(UBC_old,nscbc_ubc,ncons);
    }
    //... Stage 1 ..........................................................
    FillPatch(*this, Stemp, nghost, time, State_Type, 0, ncons);

    if (use_nscbc) {
      compute_nscbc_face_rhs(Stemp, nscbc_order == 2);
      copy_face_data( UBC_rhs1, nscbc_rhs_bc, ncons);
    }
    compute_rhs(Stemp, Real(0.5) * dt, fr_as_crse, fr_as_fine);        // Stemp = RHS(t_n,y_n) 

    MultiFab::LinComb(S2, Real(1.0), S1, 0,dt, Stemp, 0, 0, ncons, 0); // U* = U_n + dt *RHS(t_n,y_n)

    if (use_nscbc) {
      // UBC* = UBC^n + dt RHS_BC^n
      saxpy_face_data(nscbc_ubc,dt,UBC_rhs1,ncons);
    }

    //... Stage 2 ..........................................................

    // After FillPatch, Stemp = U* = U^n + dt*dUdt^n
    state[0].setNewTimeLevel(time + dt);
    FillPatch(*this, Stemp, nghost, time + dt, State_Type, 0, ncons);

    if (use_nscbc) {
      /*  nscbc_ubc now contains UBC*, so this evaluates RHS_BC(UBC*, U*).*/
      compute_nscbc_face_rhs(Stemp,nscbc_order == 2);
      copy_face_data(UBC_rhs2, nscbc_rhs_bc, ncons);
    }

    compute_rhs(Stemp, Real(0.5) * dt, fr_as_crse, fr_as_fine); // Stemp = RHS(t_n+1/2,y_*) 

    // U^(n+1/2) = 0.5*(U_n+ U*)  intermediate state
    MultiFab::LinComb(S2, Real(0.5), S1, 0,Real(0.5), S2, 0,0, ncons, 0);  

    // U^(n+1) = U^(n+1/2)  + 0.5 dt *RHS(t_n+1/2,y_*) 
    MultiFab::Saxpy(S2, Real(0.5) * dt, Stemp, 0, 0, ncons, 0);

    // NSBC BC
    if (use_nscbc) {
      /* UBC^(n+1/2) = 0.5 UBC^n + 0.5 UBC*/
      lincomb_face_data( nscbc_ubc, Real(0.5), UBC_old, Real(0.5), nscbc_ubc, ncons);
      /* UBC^(n+1)   = 0.5 UBC^n + 0.5 UBC + 0.5 dt RHS_BC(UBC*,U*) */
      saxpy_face_data( nscbc_ubc, Real(0.5)*dt, UBC_rhs2, ncons);
      update_nscbc_face_primitives();
    }
  ////////////////////////////////////////////////////////////////////////////
  } else if (order_rk == 0) {  
    // returns rhs  (debugging only)
    FillPatch(*this, Stemp, nghost, time, State_Type, 0, ncons);
    compute_rhs(Stemp, dt, fr_as_crse, fr_as_fine);
    MultiFab::Copy(S2, Stemp, 0, 0, ncons, 0);
  ////////////////////////////////////////////////////////////////////////////  
  } else if (order_rk == 1) {
    // Euler Scheme // 
    amrex::Print() << "Euler Scheme" << std::endl; 

    FillPatch(*this, Stemp, nghost, time, State_Type, 0,  ncons);  // filled at t_n to evalulate F(t_n,y_n).

    if (use_nscbc) {
      compute_nscbc_face_rhs(Stemp,nscbc_order == 2);
    }

    compute_rhs(Stemp, dt, fr_as_crse, fr_as_fine);                      // Stemp = RHS(t_n,y_n)    
    MultiFab::LinComb(S2, Real(1.0), S1, 0, dt, Stemp, 0, 0, ncons, 0);  // U_n+1 = U_n + dt *RHS(t_n,y_n) 

    if (use_nscbc) {
      saxpy_face_data(nscbc_ubc,dt,nscbc_rhs_bc,ncons);
      update_nscbc_face_primitives();
    }

  ////////////////////////////////////////////////////////////////////////////    
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
  //
  // NO SAFETY NET: if any fluid cell becomes NaN/Inf/non-positive, that is
  // a numerical failure of the scheme (under-resolved shocks, wrong CFL,
  // missing positivity limiter, etc.), not something to silently patch.
  // We detect such cells and abort with a diagnostic — upstream policy.
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
    // Pass 2: Flood-fill interior solid cells + zero their momentum.
    //
    // Runs when:
    //   - ib_move = 1 (FSI): cells that were solid can become fluid as the
    //     body moves, so solid cells must carry bounded, physically plausible
    //     data that reflects something close to wall-BC conditions.
    //   - pass2_static = 1 (opt-in for static geometry): some shock–body
    //     interaction cases benefit from clean solid-cell data because
    //     WENO stencils reach across the surface and read solid values.
    //     Empirically, flood-filling helps simple geometries
    //     (2d_bvh_cpu, airfoil_static) but *destabilises* complex geometries
    //     (2DSphere, complex_geom) — the averaged post-shock state leaks
    //     into the body and the next step's WENO oscillates on the gradient.
    //     Hence opt-in, not default.
    //
    // Zeroing momentum prevents spurious velocity amplification when
    // averaging fluid neighbors from different acoustic phases.
    //
    // Sync S2 ghost cells from neighboring fabs first: the RK stages only
    // write the valid region of S2, so cross-fab ghost cells still hold
    // data from the previous step. The flood-fill reads 3×3 neighbors,
    // and at box boundaries those neighbors land in that stale ghost
    // region — without FillBoundary, Pass 2 averages current-step valid
    // cells with previous-step ghost values, silently polluting the solid
    // state and feeding garbage into the next step's WENO stencils.
    // ------------------------------------------------------------------------
    if (CNS::ib_move || CNS::pass2_static) {
    S2.FillBoundary(geom.periodicity());
    {
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
          state(i,j,k, PC::UMZ) = Real(0.0);  // always: index exists in 2D
          state(i,j,k, PC::UET) -= ke;
        });
      }
    } // inner block
    } // end Pass 2 (ib_move || pass2_static)

    // ------------------------------------------------------------------------
    // Hard NaN/Inf/non-positive check on fluid cells. If any fluid cell is
    // broken, abort with a useful diagnostic. No silent repair — if this
    // fires, the scheme itself failed and needs fixing (resolution, CFL,
    // positivity limiter, viscosity, ...).
    //
    // When cns.strict_positivity = 1, we also abort if any fluid cell's
    // density or internal energy has collapsed to near the cons2prims
    // clipping floors (smallr ~ 1e-19, ei_min ~ 2.5e-8). That catches
    // silent clipping — the scheme may have kept marching because prims
    // were clamped, but the underlying conservative state is unphysical.
    // ------------------------------------------------------------------------
    {
      const bool strict = CNS::strict_positivity;
      // Threshold: 1e6 × the hard clipping floor. Well below any physical
      // value (nominal rho ~ 1, nominal eint ~ 2e5 J/kg) yet far enough
      // above the floor that numerical noise doesn't trip it.
      // Capture the floor constants as locals so CUDA device lambdas don't
      // reach back into CNSConstants:: namespace storage.
      const Real smallr_local = CNSConstants::smallr;
      const Real rho_floor = smallr_local * Real(1.0e6);
      const Real ei_floor  = cls_h.get_ei_min() * Real(1.0e6);

      ReduceOps<ReduceOpSum, ReduceOpSum, ReduceOpMin, ReduceOpMin> rop;
      ReduceData<int, int, Real, Real> rdata(rop);
      using RT = typename decltype(rdata)::Type;

      for (MFIter mfi(S2, false); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        Array4<Real> const& state = S2.array(mfi);
        const auto& ibMarkers = ib_mf.array(mfi);

        rop.eval(bx, rdata, [=] AMREX_GPU_DEVICE(int i, int j, int k) -> RT {
          using PC = PROB::ProbClosures;
          if (ibMarkers(i,j,k,0) != 0) return {0, 0, Real(1e300), Real(1e300)};
          const Real rho = state(i,j,k, PC::URHO);
          const Real E   = state(i,j,k, PC::UET);
          const int nonfinite = (!std::isfinite(rho) || !std::isfinite(E)) ? 1 : 0;
          int nonpos = (rho <= Real(0.0) || E <= Real(0.0)) ? 1 : 0;
          if (strict && !nonfinite && !nonpos) {
            // Approaching-clipping check (uses conservative E, not eint).
            // eint ≈ E/rho - 0.5*|u|² ; comparing E/rho against ei_floor
            // is a lower bound (true eint is smaller when KE > 0, so this
            // is conservative w.r.t. abort).
            const Real eint_approx = E / amrex::max(rho, smallr_local);
            if (rho < rho_floor || eint_approx < ei_floor) nonpos = 1;
          }
          return {nonfinite, nonpos, rho, E};
        });
      }
      auto hv = rdata.value(rop);
      int  n_nonfinite = amrex::get<0>(hv);
      int  n_nonpos    = amrex::get<1>(hv);
      Real rho_min     = amrex::get<2>(hv);
      Real E_min       = amrex::get<3>(hv);
      ParallelDescriptor::ReduceIntSum(n_nonfinite);
      ParallelDescriptor::ReduceIntSum(n_nonpos);
      ParallelDescriptor::ReduceRealMin(rho_min);
      ParallelDescriptor::ReduceRealMin(E_min);

      if (n_nonfinite > 0 || n_nonpos > 0) {
        const char* badlabel = strict
          ? "rho<=0, E<=0, rho<rho_floor, or eint<ei_floor (strict_positivity=1)"
          : "rho<=0 or E<=0";
        amrex::Print() << "\n========================================================\n"
                       << "NUMERICAL FAILURE at level " << level
                       << ", t = " << time + dt << ", dt = " << dt << "\n"
                       << "  fluid cells with NaN/Inf : " << n_nonfinite << "\n"
                       << "  fluid cells with " << badlabel << " : " << n_nonpos << "\n"
                       << "  min(rho) = " << rho_min
                       << "   (rho_floor = " << rho_floor << ")\n"
                       << "  min(E)   = " << E_min
                       << "   (ei_floor*min(rho) ~ " << ei_floor * amrex::max(rho_min, Real(0.0)) << ")\n"
                       << "========================================================\n";

        // Second pass: locate bad cells and print (i,j,k, x,y,z, cons+neighbors)
        const auto dx = geom.CellSizeArray();
        const auto plo = geom.ProbLoArray();
        int printed = 0;
        const int MAX_PRINT = 8;
        for (MFIter mfi(S2, false); mfi.isValid(); ++mfi) {
          const Box& bx = mfi.tilebox();
          Array4<Real> const& state = S2.array(mfi);
          const auto& ibMarkers = ib_mf.array(mfi);
          const int lo0 = bx.smallEnd(0), lo1 = bx.smallEnd(1);
          const int hi0 = bx.bigEnd(0),   hi1 = bx.bigEnd(1);
#if (AMREX_SPACEDIM == 3)
          const int lo2 = bx.smallEnd(2), hi2 = bx.bigEnd(2);
#else
          const int lo2 = 0, hi2 = 0;
#endif
          for (int k = lo2; k <= hi2 && printed < MAX_PRINT; ++k) {
          for (int j = lo1; j <= hi1 && printed < MAX_PRINT; ++j) {
          for (int i = lo0; i <= hi0 && printed < MAX_PRINT; ++i) {
            using PC = PROB::ProbClosures;
            if (ibMarkers(i,j,k,0) != 0) continue;
            const Real rho = state(i,j,k, PC::URHO);
            const Real E   = state(i,j,k, PC::UET);
            const bool bad = !std::isfinite(rho) || !std::isfinite(E)
                           || rho <= Real(0.0)   || E   <= Real(0.0);
            if (!bad) continue;
            const Real x = plo[0] + (i + Real(0.5)) * dx[0];
            const Real y = plo[1] + (j + Real(0.5)) * dx[1];
            amrex::AllPrint() << "[BAD CELL #" << printed << "] level=" << level
                              << " (i,j,k)=(" << i << "," << j << "," << k << ")"
                              << " (x,y)=(" << x << "," << y << ")\n"
                              << "  center:  rho=" << rho
                              << "  mx=" << state(i,j,k,PC::UMX)
                              << "  my=" << state(i,j,k,PC::UMY)
                              << "  E=" << E
                              << "  ibm0=" << int(ibMarkers(i,j,k,0))
                              << "  ibm1=" << int(ibMarkers(i,j,k,1)) << "\n";
            // 3x3 neighborhood dump
            amrex::AllPrint() << "  rho 3x3 (j-1..j+1, i-1..i+1):\n";
            for (int dj = 1; dj >= -1; --dj) {
              amrex::AllPrint() << "    ";
              for (int di = -1; di <= 1; ++di) {
                const int ii = i+di, jj = j+dj;
                amrex::AllPrint() << std::setw(13) << std::setprecision(4) << state(ii,jj,k,PC::URHO)
                                  << "[" << int(ibMarkers(ii,jj,k,0)) << "] ";
              }
              amrex::AllPrint() << "\n";
            }
            printed++;
          }}}
        }
        amrex::Print() << "\nFix the root cause — do not silently patch.\n"
                       << "========================================================\n";
        amrex::Abort("advance: non-finite / non-positive fluid state detected");
      }
    }
  }
#endif  // AMREX_USE_GPIBM

  return dt;
}
