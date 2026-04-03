#ifndef EIB_DATA_H_
#define EIB_DATA_H_

// ============================================================================
// eib_data.h — Data structures, constants, and helper types for the EIB solver
//
// Contains:
//   1. eib_detail namespace : SFINAE helpers for wall-model dispatch
//   2. ipow()              : Constexpr integer power
//   3. Constants           : Dimension indices, thresholds, image-point factors
//   4. gpData_t            : SoA storage for ghost-point data
//   5. surfImp_t           : SoA storage for surface image-point data
//   6. surfPhys_t          : SoA storage for reconstructed surface fields
//   7. is_gpData_t         : Type trait detecting gpData_t
//   8. FaceCSR             : CSR structure for per-FAB face iteration
//   9. CheckMode           : Interpolation stencil check policy
// ============================================================================

#include <limits>
#include <type_traits>
#include <utility>

#include <IBMultiFab.h>
#include <AMReX_GpuContainers.H>
#include <AMReX_IntVect.H>

#ifdef AMREX_USE_CGAL
#include "eib_cgal.h"
#else
#include "eib_bvh.h"
#endif

// ============================================================================
// 1. SFINAE helpers for wall-model dispatch
// ============================================================================

namespace eib_detail {
template <class...>
using void_t = void;

template <template <class...> class Op, class, class... Args>
struct is_detected_impl : std::false_type {};

template <template <class...> class Op, class... Args>
struct is_detected_impl<Op, void_t<Op<Args...>>, Args...> : std::true_type {};

template <template <class...> class Op, class... Args>
constexpr bool is_detected_v = is_detected_impl<Op, void_t<Op<Args...>>, Args...>::value;

template <class WM, class... Args>
using compute_surfIB_expr = decltype(WM::compute_surfIB(std::declval<Args>()...));

/// nvcc workaround: dispatch compute_surfIB outside constexpr-if in __device__ lambda
template <int eorder, class wallmodel_t, class cls_type, class PointArr, class Vec1D, class Prims2D>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void dispatch_compute_surfIB(
    const PointArr& xyz,
    const Vec1D& nvec,
    const Vec1D& t1vec,
    const Vec1D& t2vec,
    Prims2D& primsNormal,
    int type_solid_bc,
    const cls_type* cls)
{
    constexpr bool has_full = is_detected_v<
        compute_surfIB_expr,
        wallmodel_t,
        const Vec1D&,
        const Vec1D&,
        const Vec1D&,
        const Vec1D&,
        Prims2D&,
        const int,
        const cls_type*>;
    if constexpr (has_full) {
        wallmodel_t::compute_surfIB(xyz, nvec, t1vec, t2vec, primsNormal, type_solid_bc, cls);
    } else {
        amrex::ignore_unused(t1vec, t2vec);
        wallmodel_t::compute_surfIB(xyz, nvec, primsNormal, type_solid_bc, cls);
    }
}

} // namespace eib_detail

// ============================================================================
// 2. Constexpr integer power
// ============================================================================

AMREX_GPU_HOST_DEVICE constexpr int ipow(int base, int exp) {
    return (exp == 0) ? 1 : base * ipow(base, exp - 1);
}

// ============================================================================
// 3. Constants
// ============================================================================

// Index dimension
static constexpr int IDIM = AMREX_SPACEDIM - 1;     

// Box extra width for ghost point search, no more than cls_t::NGHOST - 1
static constexpr int GP_BOX_EXTRA = 0;

// Maximum number of IB geometries (for GPU GpuArray captures)
static constexpr int MAX_NGEOM = 16;

// Minimum number of valid fluid points required in the interpolation stencil
//
// For the first image point, if the number of available interpolation points is
// less than INTERP_THRESHOLD, an error will be raised and the corresponding IB
// information will be reported.
//
// For all subsequent image points, if the number of interpolation points is
// insufficient, the image point will be discarded and the order of the
// extrapolation will be reduced accordingly.
// In 2D (4 points total), we require 2.
// In 3D (8 points total), we require 3.
#if (AMREX_SPACEDIM == 2)
static constexpr int INTERP_THRESHOLD_GP   = 2;
static constexpr int INTERP_THRESHOLD_SURF = 1; // for surface data reconstruction, less strict
#else
static constexpr int INTERP_THRESHOLD_GP   = 3;
static constexpr int INTERP_THRESHOLD_SURF = 2; // for surface data reconstruction, less strict
#endif

// Number of attempts for the first image point placement
static constexpr int N_ATTEMPTS_GP   = 3;  
static constexpr int N_ATTEMPTS_SURF = 5;   

// ============================================================================
// 4. gpData_t — Ghost-point data (SoA)
// ============================================================================

template <int eorder_tparm, int iorder_tparm>
struct gpData_t {
  // CPU only attributes
  gpData_t() : ngps(0) {}
  int ngps;    

  // ideal number of interpolation points for each image point
  static constexpr int  N_InterP = ipow(iorder_tparm + 1, AMREX_SPACEDIM);

  // GPU/CPU attributes
  // Ghost point data
  Gpu::ManagedVector<Array1D< int, 0, IDIM>> gp_ijk;                        // Ghost point indices
  Gpu::ManagedVector<Array1D<Real, 0, IDIM>> ib_xyz;                        // IB point coordinates
  Gpu::ManagedVector<Real> disGP;                                           // Distance from IB point to ghost point
    
  // Surface identification
  Gpu::ManagedVector<int> geomIdx;                                          // Geometry index
  Gpu::ManagedVector<int> elemIdx;                                          // face/edge element index

  // Image point data arrays
  Gpu::ManagedVector<Array2D<Real, 0, eorder_tparm - 1, 0, IDIM>> imp_xyz;  // Image point coordinates
  Gpu::ManagedVector<Array2D< int, 0, eorder_tparm - 1, 0, IDIM>> imp_ijk;  // Image point indices
  Gpu::ManagedVector<Array1D<Real, 0, eorder_tparm - 1>> disIM;             // Distance from IB point to image points
      
  // Interpolation data for image points
  Gpu::ManagedVector<Array1D< int, 0, eorder_tparm - 1>> imp_ninterp;       // Actual number of interpolation points used for each image point
  Gpu::ManagedVector<Array3D< int, 0, eorder_tparm - 1, 0, N_InterP - 1, 0, IDIM>> imp_ip_ijk;
  Gpu::ManagedVector<Array2D<Real, 0, eorder_tparm - 1, 0, N_InterP - 1>> imp_ipweights;

  // Helper function to resize all vectors
  void resize(int n) {
      ngps = n;
      
      gp_ijk.resize(n);
      ib_xyz.resize(n);
      disGP.resize(n);
      
      geomIdx.resize(n);
      elemIdx.resize(n);

      imp_xyz.resize(n);
      imp_ijk.resize(n);
      disIM.resize(n);
      
      imp_ninterp.resize(n);
      imp_ip_ijk.resize(n);
      imp_ipweights.resize(n);
  }

  // Helper function to reserve memory for all vectors
  void reserve(int n) {
      gp_ijk.reserve(n);
      ib_xyz.reserve(n);
      disGP.reserve(n);
      
      geomIdx.reserve(n);
      elemIdx.reserve(n);

      imp_xyz.reserve(n);
      imp_ijk.reserve(n);
      disIM.reserve(n);
      
      imp_ninterp.reserve(n);
      imp_ip_ijk.reserve(n);
      imp_ipweights.reserve(n);
  }

  // Clear and free memory
  void clear() {
      ngps = 0;

      gp_ijk.clear();      
      ib_xyz.clear();      
      disGP.clear();       

      geomIdx.clear();     
      elemIdx.clear();

      imp_xyz.clear();     
      imp_ijk.clear();     
      disIM.clear();       

      imp_ninterp.clear(); 
      imp_ip_ijk.clear();  
      imp_ipweights.clear(); 
  }

  // Explicitly release memory
  void shrink() {
      gp_ijk.shrink_to_fit();
      ib_xyz.shrink_to_fit();
      disGP.shrink_to_fit();
      
      geomIdx.shrink_to_fit();
      elemIdx.shrink_to_fit();

      imp_xyz.shrink_to_fit();
      imp_ijk.shrink_to_fit();
      disIM.shrink_to_fit();
      
      imp_ninterp.shrink_to_fit();
      imp_ip_ijk.shrink_to_fit();
      imp_ipweights.shrink_to_fit();
  }
};

// ============================================================================
// 5. surfImp_t — Surface image-point data (SoA)
// ============================================================================

template <int eorder_tparm_surf, int iorder_tparm_surf>
struct surfImp_t {
  // ideal number of interpolation points for each image point
  static constexpr int  N_InterP = ipow(iorder_tparm_surf + 1, AMREX_SPACEDIM);

  // Surface identification
  Gpu::ManagedVector<int> elemIdx;      // Global face index across all geometries

  // Image point data (per face)
  Gpu::ManagedVector<Array2D<Real, 0, eorder_tparm_surf - 1, 0, IDIM>> imp_xyz;                       // Physical-space coordinates of image points placed along the outward normal
  Gpu::ManagedVector<Array2D< int, 0, eorder_tparm_surf - 1, 0, IDIM>> imp_ijk;                       // Index of the "bottom-left" grid cell associated with each image point
  Gpu::ManagedVector<Array1D<Real, 0, eorder_tparm_surf - 1>> disIM;                                  // Normal distances from the surface (IB point) to each image point
  
  // Interpolation data for image points (per face)
  Gpu::ManagedVector<Array1D< int, 0, eorder_tparm_surf - 1>> imp_ninterp;                            // Actual number of interpolation points used for each image point
  Gpu::ManagedVector<Array3D< int, 0, eorder_tparm_surf - 1, 0, N_InterP - 1, 0, IDIM>> imp_ip_ijk;   // Indices of the 8-point interpolation stencil for each image point
  Gpu::ManagedVector<Array2D<Real, 0, eorder_tparm_surf - 1, 0, N_InterP - 1>> imp_ipweights;         // Trilinear interpolation weights for the 8-point stencil of each image point

  void resize(int n) {
      elemIdx.resize(n);
      imp_xyz.resize(n);
      imp_ijk.resize(n);
      disIM.resize(n);
      imp_ninterp.resize(n);
      imp_ip_ijk.resize(n);
      imp_ipweights.resize(n);
  }

  void clear() {
      elemIdx.clear();
      imp_xyz.clear();
      imp_ijk.clear();
      disIM.clear();
      imp_ninterp.clear();
      imp_ip_ijk.clear();
      imp_ipweights.clear();
  }

  void shrink() {
      // Only shrink if capacity is significantly larger than size (e.g. > 4x)
      // to avoid frequent reallocations (memory jitter).
      if (elemIdx.capacity() > static_cast<std::size_t>(4 * elemIdx.size())) {
          elemIdx.shrink_to_fit();
          imp_xyz.shrink_to_fit();
          imp_ijk.shrink_to_fit();
          disIM.shrink_to_fit();
          imp_ninterp.shrink_to_fit();
          imp_ip_ijk.shrink_to_fit();
          imp_ipweights.shrink_to_fit();
      }
  }
};

// ============================================================================
// 6. surfPhys_t — Surface physical data (SoA)
// ============================================================================

struct surfPhys_t {

  // CPU only attributes
  surfPhys_t() : filled_elems(0) {}
  int filled_elems;

  // Indexing info
  Gpu::ManagedVector<int> elemIdx;      // Global face index across all geometries
  Gpu::ManagedVector<int> ifab;         // Local FAB index on this MPI rank
  Gpu::ManagedVector<int> lev;          // AMR level
  Gpu::ManagedVector<int> rank;         // Owning MPI rank
  Gpu::ManagedVector<int> elemfound;    // Whether this face has been located (int for GPU compatibility)
  Gpu::ManagedVector<int> ip_quality;   // number of fluid interpolation points used for the first image point
  //elemfound: 0 = outside this level, 1 = owned/valid, -1 = reserved by other rank

  // Surface fields (per face)
  Gpu::ManagedVector<Real> pressure;    // reconstructed local surface pressure
  Gpu::ManagedVector<Real> tau1;        // reconstructed local surface shear stress 1
  Gpu::ManagedVector<Real> tau2;        // reconstructed local surface shear stress 2
  Gpu::ManagedVector<Real> temperature; // reconstructed temperature
  Gpu::ManagedVector<Real> dTdn;        // reconstructed grad(T)·n
  
  // Helper function to resize all vectors
  // Note: resize() initializes new elements to 0 / default constructor.
  // If you need specific default values (e.g. -1 for indices), set them manually after resize.
  void resize(int n) {
    int old_n = elemIdx.size();
    
    elemIdx.resize(n);
    ifab.resize(n);
    lev.resize(n);
    rank.resize(n);
    elemfound.resize(n);
    
    pressure.resize(n);
    tau1.resize(n);
    tau2.resize(n);
    temperature.resize(n);
    dTdn.resize(n);
    ip_quality.resize(n);

    // Initialize new elements with specific defaults if n > old_n
    if (n > old_n) {
        for (int i = old_n; i < n; ++i) {
            ifab[i] = -1;
            lev[i]  = -1;
            rank[i] = -99;
            elemfound[i] = 0; // false
            ip_quality[i] = -1;
        }
    }
  }

  // Clear and free memory
  void clear() {
      filled_elems = 0;
      
      elemIdx.clear();
      ifab.clear();    
      lev.clear();     
      rank.clear();   
      elemfound.clear();
      
      pressure.clear(); 
      tau1.clear();     
      tau2.clear();     
      temperature.clear(); 
      dTdn.clear();     
      ip_quality.clear();     
  }

  // Explicitly release memory
  void shrink() {
    if (elemIdx.capacity() > static_cast<std::size_t>(4 * elemIdx.size())) {
        elemIdx.shrink_to_fit();
        ifab.shrink_to_fit();
        lev.shrink_to_fit();
        rank.shrink_to_fit();
        elemfound.shrink_to_fit();
        
        pressure.shrink_to_fit();
        tau1.shrink_to_fit();
        tau2.shrink_to_fit();
        temperature.shrink_to_fit();
        dTdn.shrink_to_fit();
        ip_quality.shrink_to_fit();
    }
  }

  // Reset metadata for regrid
  void reset() {
    int n = elemIdx.size();
    for (int i = 0; i < n; ++i) {
        ifab[i] = -1;
        lev[i]  = -1;
        rank[i] = -1;
        elemfound[i] = 0; // false
        ip_quality[i] = -1;
    }
  }

  AMREX_FORCE_INLINE
  bool owned(int f_idx, int lev, int ifab) const noexcept
  {
        return elemfound[f_idx] == 1 &&
                this->lev[f_idx]  == lev &&
                this->ifab[f_idx] == ifab;
  }

};

// ============================================================================
// 7. Type traits and helpers
// ============================================================================

// Type trait to detect if a type is gpData_t (has gp_ijk member)
template <typename T, typename = void>
struct is_gpData_t : std::false_type {};

template <typename T>
struct is_gpData_t<T, std::void_t<decltype(std::declval<T>().gp_ijk)>> : std::true_type {};

// ============================================================================
// 8. GPStoreView / GPStore — Level-wide flattened ghost-point storage (CSR)
// ============================================================================

/// \brief GPU-capturable POD view into GPStore.  Holds raw pointers only.
///        This struct is trivially copyable and can be captured by GPU lambdas.
template <int eorder_tparm, int iorder_tparm>
struct GPStoreView {
  static constexpr int N_InterP = ipow(iorder_tparm + 1, AMREX_SPACEDIM);

  int  total_ngps;                                   // total ghost points on this level
  int  nfabs;                                        // number of local FABs
  const int*  fab_offsets;                            // CSR offsets: fab_offsets[ifab] .. fab_offsets[ifab+1]

  // Per-GP arrays (indexed 0..total_ngps-1)
  const int*                    gp_fab;                   // local FAB index for each GP
  const Array1D< int, 0, IDIM>* gp_ijk;
  const Array1D<Real, 0, IDIM>* ib_xyz;
  const Real*                   disGP;
  const int*                    geomIdx;
  const int*                    elemIdx;
  const Array2D<Real, 0, eorder_tparm - 1, 0, IDIM>* imp_xyz;
  const Array2D< int, 0, eorder_tparm - 1, 0, IDIM>* imp_ijk;
  const Array1D<Real, 0, eorder_tparm - 1>*           disIM;
  const Array1D< int, 0, eorder_tparm - 1>*           imp_ninterp;
  const Array3D< int, 0, eorder_tparm - 1, 0, N_InterP - 1, 0, IDIM>* imp_ip_ijk;
  const Array2D<Real, 0, eorder_tparm - 1, 0, N_InterP - 1>*          imp_ipweights;

  // Non-const data pointers for initialiseGPs (write pass)
  Array1D< int, 0, IDIM>* gp_ijk_w;
  Array1D<Real, 0, IDIM>* ib_xyz_w;
  Real*                   disGP_w;
  int*                    geomIdx_w;
  int*                    elemIdx_w;
  Array2D<Real, 0, eorder_tparm - 1, 0, IDIM>* imp_xyz_w;
  Array2D< int, 0, eorder_tparm - 1, 0, IDIM>* imp_ijk_w;
  Array1D<Real, 0, eorder_tparm - 1>*           disIM_w;
  Array1D< int, 0, eorder_tparm - 1>*           imp_ninterp_w;
  Array3D< int, 0, eorder_tparm - 1, 0, N_InterP - 1, 0, IDIM>* imp_ip_ijk_w;
  Array2D<Real, 0, eorder_tparm - 1, 0, N_InterP - 1>*          imp_ipweights_w;

  /// Get GP range for a local FAB index
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  int gp_begin(int ifab) const { return fab_offsets[ifab]; }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  int gp_end(int ifab) const { return fab_offsets[ifab + 1]; }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  int ngps(int ifab) const { return fab_offsets[ifab + 1] - fab_offsets[ifab]; }
};

/// \brief Level-wide flattened ghost-point storage.  Owns memory via ManagedVector.
///        CSR layout: fab_offsets[ifab] gives the first GP index for local fab ifab.
template <int eorder_tparm, int iorder_tparm>
struct GPStore {
  static constexpr int N_InterP = ipow(iorder_tparm + 1, AMREX_SPACEDIM);

  int total_ngps = 0;   // total ghost points on this level (sum over all FABs)
  int nfabs      = 0;   // number of local FABs

  // CSR index: size = nfabs + 1.  fab_offsets[i] = start index of FAB i's GPs.
  Gpu::ManagedVector<int> fab_offsets;

  // Per-GP arrays (flattened, size = total_ngps each)
  Gpu::ManagedVector<int>                    gp_fab;   // gp_fab[ii] = local FAB index for GP ii
  Gpu::ManagedVector<Array1D< int, 0, IDIM>> gp_ijk;
  Gpu::ManagedVector<Array1D<Real, 0, IDIM>> ib_xyz;
  Gpu::ManagedVector<Real>                   disGP;
  Gpu::ManagedVector<int>                    geomIdx;
  Gpu::ManagedVector<int>                    elemIdx;
  Gpu::ManagedVector<Array2D<Real, 0, eorder_tparm - 1, 0, IDIM>> imp_xyz;
  Gpu::ManagedVector<Array2D< int, 0, eorder_tparm - 1, 0, IDIM>> imp_ijk;
  Gpu::ManagedVector<Array1D<Real, 0, eorder_tparm - 1>>           disIM;
  Gpu::ManagedVector<Array1D< int, 0, eorder_tparm - 1>>           imp_ninterp;
  Gpu::ManagedVector<Array3D< int, 0, eorder_tparm - 1, 0, N_InterP - 1, 0, IDIM>> imp_ip_ijk;
  Gpu::ManagedVector<Array2D<Real, 0, eorder_tparm - 1, 0, N_InterP - 1>>          imp_ipweights;

  /// Allocate flat arrays from per-fab counts.
  /// \param counts  Vector of GP counts per local FAB (size = nfabs_in).
  void allocate(const amrex::Vector<int>& counts) {
    nfabs = static_cast<int>(counts.size());

    // Build CSR offsets via exclusive prefix sum
    fab_offsets.resize(nfabs + 1);
    fab_offsets[0] = 0;
    for (int i = 0; i < nfabs; ++i) {
      fab_offsets[i + 1] = fab_offsets[i] + counts[i];
    }
    total_ngps = fab_offsets[nfabs];

    // Resize all per-GP arrays
    gp_fab.resize(total_ngps);
    gp_ijk.resize(total_ngps);

    // Fill gp_fab: expand CSR offsets to per-GP FAB index
    for (int f = 0; f < nfabs; ++f) {
      for (int g = fab_offsets[f]; g < fab_offsets[f + 1]; ++g) {
        gp_fab[g] = f;
      }
    }
    ib_xyz.resize(total_ngps);
    disGP.resize(total_ngps);
    geomIdx.resize(total_ngps);
    elemIdx.resize(total_ngps);
    imp_xyz.resize(total_ngps);
    imp_ijk.resize(total_ngps);
    disIM.resize(total_ngps);
    imp_ninterp.resize(total_ngps);
    imp_ip_ijk.resize(total_ngps);
    imp_ipweights.resize(total_ngps);
  }

  /// Return a GPU-capturable read/write view of this store.
  GPStoreView<eorder_tparm, iorder_tparm> view() const {
    GPStoreView<eorder_tparm, iorder_tparm> v;
    v.total_ngps   = total_ngps;
    v.nfabs        = nfabs;
    v.fab_offsets   = fab_offsets.data();

    v.gp_fab        = gp_fab.data();
    v.gp_ijk        = gp_ijk.data();
    v.ib_xyz        = ib_xyz.data();
    v.disGP         = disGP.data();
    v.geomIdx       = geomIdx.data();
    v.elemIdx       = elemIdx.data();
    v.imp_xyz       = imp_xyz.data();
    v.imp_ijk       = imp_ijk.data();
    v.disIM         = disIM.data();
    v.imp_ninterp   = imp_ninterp.data();
    v.imp_ip_ijk    = imp_ip_ijk.data();
    v.imp_ipweights = imp_ipweights.data();

    // const_cast for writable pointers (initialiseGPs write pass)
    v.gp_ijk_w        = const_cast<Array1D< int, 0, IDIM>*>(gp_ijk.data());
    v.ib_xyz_w        = const_cast<Array1D<Real, 0, IDIM>*>(ib_xyz.data());
    v.disGP_w         = const_cast<Real*>(disGP.data());
    v.geomIdx_w       = const_cast<int*>(geomIdx.data());
    v.elemIdx_w       = const_cast<int*>(elemIdx.data());
    v.imp_xyz_w       = const_cast<Array2D<Real, 0, eorder_tparm - 1, 0, IDIM>*>(imp_xyz.data());
    v.imp_ijk_w       = const_cast<Array2D< int, 0, eorder_tparm - 1, 0, IDIM>*>(imp_ijk.data());
    v.disIM_w         = const_cast<Array1D<Real, 0, eorder_tparm - 1>*>(disIM.data());
    v.imp_ninterp_w   = const_cast<Array1D< int, 0, eorder_tparm - 1>*>(imp_ninterp.data());
    v.imp_ip_ijk_w    = const_cast<Array3D< int, 0, eorder_tparm - 1, 0, N_InterP - 1, 0, IDIM>*>(imp_ip_ijk.data());
    v.imp_ipweights_w = const_cast<Array2D<Real, 0, eorder_tparm - 1, 0, N_InterP - 1>*>(imp_ipweights.data());

    return v;
  }

  void clear() {
    total_ngps = 0;
    nfabs      = 0;
    fab_offsets.clear();
    gp_fab.clear();
    gp_ijk.clear();
    ib_xyz.clear();
    disGP.clear();
    geomIdx.clear();
    elemIdx.clear();
    imp_xyz.clear();
    imp_ijk.clear();
    disIM.clear();
    imp_ninterp.clear();
    imp_ip_ijk.clear();
    imp_ipweights.clear();
  }

  void shrink() {
    auto cap = gp_ijk.capacity();
    if (cap > static_cast<std::size_t>(4 * total_ngps) && cap > 1000u) {
      fab_offsets.shrink_to_fit();
      gp_fab.shrink_to_fit();
      gp_ijk.shrink_to_fit();
      ib_xyz.shrink_to_fit();
      disGP.shrink_to_fit();
      geomIdx.shrink_to_fit();
      elemIdx.shrink_to_fit();
      imp_xyz.shrink_to_fit();
      imp_ijk.shrink_to_fit();
      disIM.shrink_to_fit();
      imp_ninterp.shrink_to_fit();
      imp_ip_ijk.shrink_to_fit();
      imp_ipweights.shrink_to_fit();
    }
  }
};

// ============================================================================
// 9. FaceCSR — CSR structure for per-FAB face iteration
// ============================================================================

struct FaceCSR {
  Gpu::ManagedVector<int> fab_offsets;   // Offsets for each FAB (size = nfab + 1)
  Gpu::ManagedVector<int> face_indices;  // Contiguous array of face indices

  void clear() {
        fab_offsets.clear();  
        face_indices.clear();
  }

  void shrink() {
        fab_offsets.shrink_to_fit();
        face_indices.shrink_to_fit();
  }
};

// ============================================================================
// 10. CheckMode — Interpolation stencil check policy
// ============================================================================

enum class CheckMode {
    Silent,      // Do not output anything, just return status
    Warn,        // Output a warning message, return status
    Abort        // Abort execution immediately on failure
};

#endif // EIB_DATA_H_
