#ifndef EBMultiFab_H_
#define EBMultiFab_H_

#include <AMReX_Amr.H>
#include <AMReX_AmrLevel.H>
#include <AMReX_FabArray.H>

///
/// \brief EBFab holds the solid point and ghost point boolean arrays
/// \param marker_t type of marker field (bool or Real) 
///
template<typename marker_t>
class EBFab: public BaseFab<bool> {
 public:

  // using Box
  explicit inline EBFab<marker_t>(const Box& b, int ncomp, bool alloc = true, bool shared = false, Arena* ar = nullptr) : BaseFab<marker_t>(b, ncomp, alloc, shared, ar) {};
  // using EBFab
  explicit inline EBFab<marker_t>(const EBFab<marker_t>& rhs, MakeType make_type, int scomp, int ncomp) : BaseFab<marker_t>(rhs, make_type, scomp, ncomp) {};

  ~EBFab() {};
};
///
/// \brief EBMultiFab holds an array of IBFab on a level 
///
template<typename marker_t>
class EBMultiFab : public FabArray<EBFab<marker_t>> {
 public:
  // constructor from BoxArray and DistributionMapping
  explicit EBMultiFab(
      const BoxArray& bxs, const DistributionMapping& dm, const int nvar,
      const int ngrow, const MFInfo& info = MFInfo{true, amrex::The_Managed_Arena()},
      const FabFactory<EBFab<marker_t>>& factory = DefaultFabFactory<EBFab<marker_t>>()) : FabArray<EBFab<marker_t>>(bxs, dm, nvar, ngrow, info, factory) {};

  EBMultiFab(EBMultiFab&& rhs) noexcept 
                        : FabArray<EBFab<marker_t>>(std::move(rhs)) {};

  ~EBMultiFab() {};

  /////////////////////////////////////////////////////////////////
  void inline copytoRealMF(MultiFab& mf, int ebcomp, int mfcomp) {
    for (MFIter mfi(*this, false); mfi.isValid(); ++mfi) {      
      const Box& ebbox = mfi.validbox();  // box without ghost points

      EBFab<marker_t>& ebfab = this->get(mfi);      
      const Array4<bool>& ebMarkers = ebfab.array();    // boolean array
      const Array4<Real>& realfield = mf.array(mfi);  // real array

      amrex::ParallelFor(
          ebbox, 1, [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
            realfield(i, j, k, mfcomp + n) = ebMarkers(i, j, k, ebcomp + n);
          });
    }
  }

};

#endif