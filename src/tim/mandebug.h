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
// e.g : checkNaN_prims(" check flux array",mfi,cls_h.NCONS,flx_x);    
void inline checkNaN_prims(const std::string& errorMessage, const MFIter& mfi,const int NVAR,
                           const Array4<Real>& prims)
{

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
//------------------------------------------------------------------//
// print info at  a point
// e.g. use: 
// printinfo_point(" point check 11 20",mfi, 11,20,0,cls_h.NCONS,cls_h.NPRIM, prims,cons,rhs);   
void inline printinfo_point(const std::string& errorMessage,const MFIter& mfi,
  const int i,const int j , const int k, const int NCONS, const int NPRIM,
  const Array4<Real>& prims,const Array4<Real>& cons,const Array4<Real>& rhs)
{
  const Box& box  = mfi.growntilebox(0);
  amrex::ParallelFor(
  box, [=] AMREX_GPU_DEVICE(int is, int js, int ks) noexcept
  {    
      if ((i==is) && (j==js) & (k==ks)) {
        std::cout << "  ------------------------------------- \n";
        std::cout << " printvalue at :"  <<  errorMessage <<  "\n";
        printf(" i=%d  j=%d k=%d \n",i,j,k);
        for (int n = 0; n < NCONS; n++) {   
          printf(" n=%d U=%f RHS=%f \n",n,cons(i,j,k,n),rhs(i,j,k,n));
        }
        printf(" -- \n");
        for (int n = 0; n < NPRIM; n++) {   
          printf(" n=%d Q=%f \n",n,prims(i,j,k,n));
        }
        std::cout << "  ------------------------------------- \n";
      }
  });  
}  
//------------------------------------------------------------------//
// print array at a known ijk point
// e.g. use: printarray_point(" show point ",i.j,k,cls_t::NPRIMS,prims);
void inline printarray_point(const std::string& errorMessage,
  const int i,const int j , const int k, const int NPRIM,const Array4<Real>& prims)
{
  std::cout << "  ------------------------------------- \n";
  std::cout << " printvalue at :"  <<  errorMessage <<  "\n";
  printf(" i=%d  j=%d k=%d \n",i,j,k);
  printf(" -- \n");
  for (int n = 0; n < NPRIM; n++) {   
    printf(" n=%d Q=%f \n",n,prims(i,j,k,n));
  }
  std::cout << "  ------------------------------------- \n";    
}  
//------------------------------------------------------------------//
// print flux info at  a point
//  e.g:   print_flx_point(" in flux X ",mfi,37,71,0, cls_h.NCONS,flx,cons,state);
void inline print_flx_point(const std::string& errorMessage,const MFIter& mfi,
  const int i,const int j , const int k, const int NCONS,
  const Array4<Real>& flx,const Array4<Real>& cons,const Array4<Real>& rhs)
{
  const Box& box  = mfi.growntilebox(0);
  amrex::ParallelFor(
  box, [=] AMREX_GPU_DEVICE(int is, int js, int ks) noexcept
  {    
      if ((i==is) && (j==js) & (k==ks)) {
        std::cout << "  ------------------------------------- \n";
        std::cout << " printfluxvalue at :"  <<  errorMessage <<  "\n";
        printf(" i=%d  j=%d k=%d \n",i,j,k);
        for (int n = 0; n < NCONS; n++) {   
          printf(" n=%d U=%f RHS=%f \n",n,cons(i,j,k,n),rhs(i,j,k,n));
        }
        printf(" -- \n");
        for (int n = 0; n < NCONS; n++) {   
          printf(" n=%d F(i)=%f F(i+1)=%f \n",n,flx(i,j,k,n),flx(i+1,j,k,n));
        }
        std::cout << "  ------------------------------------- \n";
      }
  });  
}  

#endif




