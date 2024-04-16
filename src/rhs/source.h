#ifndef SOURCE_H_
#define SOURCE_H_

#include <prob.h>

using namespace amrex;
//using namespace PROB;

// manual source
template <typename cls_t>
class user_source_t {
  public:
  void inline src(const amrex::MFIter &mfi,
                  const amrex::Array4<const amrex::Real> &prims,
                  const amrex::Array4<amrex::Real> &rhs, const cls_t *cls_d,
                  amrex::Real dt){

   // std:: cout << " SNM: dt " << dt << std::endl; 

  //  Abort("src to do");

    BL_PROFILE("user_source_t::src()");

    // put here because this is a .h file
    // using amrex::Array4;
    // using amrex::Box;
    // using amrex::FArrayBox;
    // using amrex::IArrayBox;
    // using amrex::Real;

    const Box bx = mfi.tilebox();

 //   PROB::ProbParm const* lprobparm;

    //const Real *dx = geomdata.CellSize();

      amrex::ParallelFor(bx,
      [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      { 
      const auto& cls = *cls_d;

    /// me dice 'PROB' has not been declared ??  (try another*.h doesn't recognised prob.h )
    //  user_source(i,j,k,prims,rhs,cls); 

    // remove source.h from rhs??


      Real rho = prims(i, j, k, cls.QRHO);
     rhs(i,j,k,cls.UMY) -= rho; // momentum (g=1) test
      rhs(i,j,k,cls.UET) -= rho*prims(i, j, k, cls.QV);  // energy

      });

  };
};
#endif
