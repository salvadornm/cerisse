#ifndef mandebug_H_
#define mandebug_H_

#include <AMReX_FluxRegister.H>
#include <AMReX_FArrayBox.H>
#include <CNS.h>

#include <cmath>

using namespace amrex;

// set of auxiliary arrays for manual debugging
//------------------------------------------------------------------//
// check array for NaNs (NVAR is size array)
void inline checkNaN_prims(const std::string& errorMessage, const MFIter& mfi,const int NVAR,
                           const Array4<Real>& prims)
{

  std::ostream& ss = std::cout; 
  const Box& box  = mfi.growntilebox(0);

  amrex::ParallelFor( box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    for (int n = 0; n < NVAR; n++)
    {
      const Real aux = prims(i,j,k,n);
      if (amrex::isnan(aux) || (aux > 1e10 ) )
      {
        std::cout << "  ------------------------------------- \n";
        std::cout << " NaN detected  in: "  <<  errorMessage <<  "\n";
        std::cout << " i j k " << i << "," << j << "," << k << "(n=" << n << "): " << "\n";        
        std::cout << "   array:  \n";
        for (int m = 0; m < NVAR; m++) {
          std::cout << "( " << m << " )" << " = " << prims(i,j,k,m) << " \n";
        }
        std::cout << "  ------------------------------------- \n";
        amrex::Abort();
      }
    }      
  }
  );

}
//------------------------------------------------------------------//
// check  NaNs in a point
void inline checkNaN_point(const std::string& errorMessage, 
  const int i,const int j , const int k, const int NVAR, const Array4<Real>& arr)
{
  for (int n = 0; n < NVAR; n++)
  {
    const Real aux = arr(i,j,k,n);
    if (amrex::isnan(aux) || (aux > 1e10 ) )
    {
      std::cout << "  ------------------------------------- \n";
      std::cout << " NaN detected  at  point in:"  <<  errorMessage <<  "\n";
      std::cout << " i j k " << i << "," << j << "," << k << "(n=" << n << "): " << "\n";        
      std::cout << "   array:  \n";
      for (int m = 0; m < NVAR; m++) {
        std::cout << "( " << m << " )" << " = " << arr(i,j,k,m) << " \n";
      }
      std::cout << "  ------------------------------------- \n";
      amrex::Abort();
    }
  }     
}
#endif




