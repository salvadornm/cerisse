#ifndef IBMultiFab_H_
#define IBMultiFab_H_

#include <AMReX_Amr.H>
#include <AMReX_AmrLevel.H>
#include <AMReX_FabArray.H>

///
/// \brief IBFab holds the solid point and ghost point boolean arrays
/// \param marker_t type of marker field (bool or Real) 
///
template<typename marker_t, typename gp_t>
class IBFab: public BaseFab<bool> {
 public:
  gp_t gpData;  

  // using Box
  explicit inline IBFab<marker_t, gp_t>(const Box& b, int ncomp, bool alloc = true, bool shared = false, Arena* ar = nullptr) : BaseFab<marker_t>(b, ncomp, alloc, shared, ar) {};
  // using IBFab
  explicit inline IBFab<marker_t, gp_t>(const IBFab<marker_t, gp_t>& rhs, MakeType make_type, int scomp, int ncomp) : BaseFab<marker_t>(rhs, make_type, scomp, ncomp) {};

  // IBFab (IBFab&& rhs) noexcept = default;
  // IBFab& operator= (IBFab&&) noexcept = default;

  // IBFab (const IBM::IBFab&) {};
  // IBFab& operator= (const IBFab&) = delete;

  ~IBFab() {};
};


///
/// \brief IBMultiFab holds an array of IBFab on a level
/// \param NIMPS Number of image points (integer) 
///
template<typename marker_t, typename gp_t>
class IBMultiFab : public FabArray<IBFab<marker_t,gp_t>> {
 public:
  // constructor from BoxArray and DistributionMapping
  explicit inline IBMultiFab<marker_t,gp_t>(
      const BoxArray& bxs, const DistributionMapping& dm, const int nvar,
      const int ngrow, const MFInfo& info = MFInfo(),
      const FabFactory<IBFab<marker_t,gp_t>>& factory = DefaultFabFactory<IBFab<marker_t,gp_t>>()) : FabArray<IBFab<marker_t,gp_t>>(bxs, dm, nvar, ngrow, info, factory) {};

  IBMultiFab<marker_t,gp_t>(IBMultiFab<marker_t,gp_t>&& rhs) noexcept 
                        : FabArray<IBFab<marker_t,gp_t>>(std::move(rhs)) {};

  ~IBMultiFab() {};

  void inline copytoRealMF(MultiFab& mf, int ibcomp, int mfcomp) {
    for (MFIter mfi(*this, false); mfi.isValid(); ++mfi) {
      // const Box& ibbox = mfi.fabbox(); // box with ghost points
      const Box& ibbox = mfi.validbox();  // box without ghost points

      IBFab<marker_t,gp_t>& ibfab = this->get(mfi);
      // const FArrayBox& realfab = mf.get(mfi);
      const Array4<bool>& ibMarkers = ibfab.array();    // boolean array
      const Array4<Real>& realfield = mf.array(mfi);  // real array

      amrex::ParallelFor(
          ibbox, 2, [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
            realfield(i, j, k, mfcomp + n) = ibMarkers(i, j, k, ibcomp + n);
          });
    }
  }
};

#endif