#ifndef IBMultiFab_H_
#define IBMultiFab_H_

#include <algorithm>
#include <AMReX_FabArray.H>
#include <AMReX_MultiFab.H>

// ============================================================================
// IBFab: extends BaseFab<marker_t> with per-FAB ghost-point data (gpData).
//
//   BaseFab<uint8_t>   ← marker array (solid/fluid/ghost flags), GPU-managed
//   gpData_t           ← ghost-point geometry + interpolation data (ManagedVector)
//
// gpData is NOT part of BaseFab memory and does NOT participate in FillBoundary
// or MPI communication. It is computed locally per FAB by computeMarkers/initialiseGPs.
// ============================================================================

/// \brief IBFab holds marker data and ghost point data
/// \tparam marker_t Type of the marker data (typically uint8_t)
/// \tparam gp_t     Type of the ghost point data (gpData_t)
template<typename marker_t, typename gp_t>
class IBFab : public amrex::BaseFab<marker_t> {
public:
  gp_t gpData;

  // Primary ctor: allocates marker array, gpData starts empty.
  explicit IBFab(const amrex::Box& b, int ncomp,
                 bool alloc = true, bool shared = false, amrex::Arena* ar = nullptr)
      : amrex::BaseFab<marker_t>(b, ncomp, alloc, shared, ar) {}

  // MakeType ctor (alias/deep-copy): gpData is NOT carried — aliases are marker-only.
  explicit IBFab(const IBFab<marker_t, gp_t>& rhs, amrex::MakeType make_type, int scomp, int ncomp)
      : amrex::BaseFab<marker_t>(rhs, make_type, scomp, ncomp) {}

  ~IBFab() = default;

  // Prevent expensive implicit deep-copy of gpData (contains Gpu::ManagedVector arrays)
  IBFab(const IBFab&) = delete;
  IBFab& operator=(const IBFab&) = delete;

  IBFab(IBFab&&) noexcept = default;
  IBFab& operator=(IBFab&&) noexcept = default;
};

// ============================================================================
// IBMultiFab: FabArray of IBFab.  Forces The_Managed_Arena so marker data is
// accessible from both CPU and GPU.  Provides copytoRealMF() for plotfile output.
// ============================================================================

/// \brief IBMultiFab holds an array of IBFab on a level
/// \tparam marker_t Type of the marker data
/// \tparam gp_t     Type of the ghost point data
template<typename marker_t, typename gp_t>
class IBMultiFab : public amrex::FabArray<IBFab<marker_t, gp_t>> {
public:
  using Fab  = IBFab<marker_t, gp_t>;
  using Base = amrex::FabArray<Fab>;

  /// \brief Default MFInfo that routes allocation to managed (CPU+GPU) memory.
  static amrex::MFInfo DataMFInfo() {
      amrex::MFInfo info;
      info.SetArena(amrex::The_Managed_Arena());
      return info;
  }

  explicit IBMultiFab(
      const amrex::BoxArray& bxs, const amrex::DistributionMapping& dm, int nvar, int ngrow,
      const amrex::MFInfo& info = DataMFInfo(),
      const amrex::FabFactory<Fab>& factory = amrex::DefaultFabFactory<Fab>())
      : Base(bxs, dm, nvar, ngrow, info, factory) {}

  ~IBMultiFab() = default;

  IBMultiFab(IBMultiFab&&) noexcept = default;
  IBMultiFab& operator=(IBMultiFab&&) noexcept = default;

  IBMultiFab(const IBMultiFab&) = delete;
  IBMultiFab& operator=(const IBMultiFab&) = delete;

  /// \brief Copy marker components from this IBMultiFab into a Real MultiFab (e.g. for plotfile).
  /// \param mf      Destination MultiFab (Real).
  /// \param ibcomp  Starting component in this IBMultiFab.
  /// \param mfcomp  Starting component in the destination MultiFab.
  void copytoRealMF(amrex::MultiFab& mf, int ibcomp, int mfcomp) {
    const int ncomp_copy = std::min(this->nComp() - ibcomp, mf.nComp() - mfcomp);
    if (ncomp_copy <= 0) return;

    for (amrex::MFIter mfi(*this, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
      const amrex::Box& bx = mfi.tilebox();
      auto const& src = this->get(mfi).array();
      auto const& dst = mf.array(mfi);

      amrex::ParallelFor(bx, ncomp_copy,
        [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
          dst(i, j, k, mfcomp + n) = static_cast<amrex::Real>(src(i, j, k, ibcomp + n));
        });
    }
  }
};

#endif // IBMultiFab_H_