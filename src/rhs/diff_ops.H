#ifndef CNS_DIFF_OPS_H_
#define CNS_DIFF_OPS_H_

#include <AMReX_Array4.H>
#include <AMReX_REAL.H>

/**
*  @brief calculate interpolation and central derivative coefficients 2/4/6 order face centred
*  (all numbers stored here)
*/

// --------------- FACE CENTRED ----------------------  i+1/2
/// i-2    i-1    i   |  i+1  i+2   i+3
///               0   |   1                      Order 2
///          0    1   |   2    3                 Order 4 
///  0       1    2   |   3    4     5           Order 6
///  -2.5h-1.5h -0.5h | 0.5h  1.5h  2.5h
// ----------------------------------------------------

// derivative coefficients face centred  2nd order
static constexpr Real fcw2_0 =  -1.0; 
static constexpr Real fcw2_1 =   1.0; 

// derivative coefficients face centred  2nd order   (Finite Difference)
static constexpr Real fd_fcw2_0 =  -1.0; 
static constexpr Real fd_fcw2_1 =   1.0; 

// interpolation coefficients face centred 2nd order (Finite Volume)
static constexpr Real int_fcw2_0 =  0.5; 
static constexpr Real int_fcw2_1 =  0.5; 

// derivative coefficients face centred  4th order   (Finite Difference)
static constexpr Real fd_fcw4_0 =   1.0 / 24.0; 
static constexpr Real fd_fcw4_1 =  -9.0 / 8.0; 
static constexpr Real fd_fcw4_2 =   9.0 / 8.0; 
static constexpr Real fd_fcw4_3 =  -1.0 / 24.0; 

// derivative coefficients face centred  4th order   (Finite Volume)
static constexpr Real fcw4_0 =   1.0 / 12.0; 
static constexpr Real fcw4_1 =  -5.0 / 4.0; 
static constexpr Real fcw4_2 =   5.0 / 4.0; 
static constexpr Real fcw4_3 =  -1.0 / 12.0; 

// interpolation coefficients face centred 4nd order (Finite Volume)
static constexpr Real int_fcw4_0 = -1.0 / 12.0; 
static constexpr Real int_fcw4_1 =  7.0 / 12.0; 
static constexpr Real int_fcw4_2 =  7.0 / 12.0; 
static constexpr Real int_fcw4_3 = -1.0 / 12.0; 

// derivative coefficients face centred  6th order   (Finite Difference)
static constexpr Real fd_fcw6_0 =    -3.0 / 640.0; 
static constexpr Real fd_fcw6_1 =    25.0 / 384.0; 
static constexpr Real fd_fcw6_2 =   -75.0 / 64.0; 
static constexpr Real fd_fcw6_3 =    75.0 / 64.0;
static constexpr Real fd_fcw6_4 =   -25.0 / 384.0;
static constexpr Real fd_fcw6_5 =     3.0 / 640.0;

//derivative coefficients face centred  6th order   (Finite Volume) 
static constexpr Real fcw6_0 =  -1.0  / 90.0; 
static constexpr Real fcw6_1 =  +5.0  / 36.0; 
static constexpr Real fcw6_2 =  -49.0 / 36.0; 
static constexpr Real fcw6_3 =  +49.0 / 36.0; 
static constexpr Real fcw6_4 =  -5.0  / 36.0; 
static constexpr Real fcw6_5 =  +1.0  / 60.0; 

// interpolation coefficients face centred 6th order  (Finite volume) 
static constexpr Real int_fcw6_0 =  1.0  /  60.0; 
static constexpr Real int_fcw6_1 =  -8.0 / 60.0; 
static constexpr Real int_fcw6_2 =  37.0 / 60.0; 
static constexpr Real int_fcw6_3 =  37.0 / 60.0; 
static constexpr Real int_fcw6_4 =  -8.0 / 60.0; 
static constexpr Real int_fcw6_5 =  1.0  / 60.0; 

// ---------------CELL  CENTRED ----------------------
/// i-3   i-2    i-1    i     i+1  i+2   i+3
///               0     *      1                Order 2
///        0      1     *      2    3           Order 4 
///  0     1      2     *      3    4     5     Order 6
///  -3h  -2h    -h     0      h    2h    3h
// ----------------------------------------------------
// coefficients checked and correct May 2025 //

// store Order 2 coefficients cell centred
static constexpr Real ccw2_0 =  -0.5; 
static constexpr Real ccw2_1 =   0.5; 

// store Order 4 coefficients cell centred
static constexpr Real ccw4_0 =   1.0 / 12.0; 
static constexpr Real ccw4_1 =  -2.0 / 3.0; 
static constexpr Real ccw4_2 =   2.0 / 3.0; 
static constexpr Real ccw4_3 =  -1.0 / 12.0; 

// store Order 6 coefficients cell centred 
static constexpr Real ccw6_0 =  -1.0 / 60.0; 
static constexpr Real ccw6_1 =   3.0 / 20.0; 
static constexpr Real ccw6_2 =  -3.0 / 4.0; 
static constexpr Real ccw6_3 =   3.0 / 4.0; 
static constexpr Real ccw6_4 =  -3.0 / 20.0; 
static constexpr Real ccw6_5 =  +1.0 / 60.0; 

// coefficients in cell faces 
template <int order>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void calc_CDcoeffs(  
        amrex::Array1D<amrex::Real,0,order>  &INTcoef,  amrex::Array1D<amrex::Real,0,order>  &CDcoef) {

    switch (order)
    {
    case 2:
      // interpolation
      INTcoef(0) =  int_fcw2_0;
      INTcoef(1) =  int_fcw2_1;
      // derivatives coefficients
      CDcoef(0)  =  fcw2_0;
      CDcoef(1)  =  fcw2_1;      
      break;
    case 4:
      // interpolation 
      INTcoef(0) =  int_fcw4_0;  
      INTcoef(1) =  int_fcw4_1;  
      INTcoef(2) =  int_fcw4_2;  
      INTcoef(3) =  int_fcw4_3;  
      // derivatives coefficients
      CDcoef(0)  =  fcw4_0;
      CDcoef(1)  =  fcw4_1;      
      CDcoef(2)  =  fcw4_2;      
      CDcoef(3)  =  fcw4_3;      
      break;
    case 6: 
      // interpolation 
      INTcoef(0) =  int_fcw6_0;
      INTcoef(1) =  int_fcw6_1;
      INTcoef(2) =  int_fcw6_2;
      INTcoef(3) =  int_fcw6_3;
      INTcoef(4) =  int_fcw6_4;
      INTcoef(5) =  int_fcw6_5;
      // derivatives coefficients 
      CDcoef(0)   =  fcw6_0;
      CDcoef(1)   =  fcw6_1;   
      CDcoef(2)   =  fcw6_2;      
      CDcoef(3)   =  fcw6_3;
      CDcoef(4)   =  fcw6_4;
      CDcoef(5)   =  fcw6_5;              
      break;
    default: 
      amrex::Abort("Diffusive order available 2/4/6 MUST specify one of those ");      
      break;
    }
    
    return;
  }


// Normal derivative (FACE CENTRED) with variable order using central differences ////
/// ********************************************************************** ////
template <int order>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real normal_diff(
  amrex::IntVect iv, int idir, int comp, amrex::Array4<amrex::Real const> const& q,
  amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& dxinv) noexcept;

/// @brief 2nd order central difference face centred
template <>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real normal_diff<2>(
  amrex::IntVect iv, int idir, int comp, amrex::Array4<amrex::Real const> const& q,
  amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& dxinv) noexcept
{
  const auto ivn = amrex::IntVect::TheDimensionVector(idir);
  return (q(iv, comp) - q(iv - ivn, comp)) * dxinv[idir];
}

/// @brief 4th order central difference face centred
template <>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real normal_diff<4>(
  amrex::IntVect iv, int idir, int comp, amrex::Array4<amrex::Real const> const& q,
  amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& dxinv) noexcept
{
  using amrex::Real;
  const auto ivn = amrex::IntVect::TheDimensionVector(idir);

  return (fcw4_3 * q(iv + ivn, comp) + fcw4_2 * q(iv, comp) + fcw4_1 * q(iv - ivn, comp) +
          fcw4_0 * q(iv - 2 * ivn, comp)) * dxinv[idir];
}

/// @brief 6th order central difference face centred
template <>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real normal_diff<6>(
  amrex::IntVect iv, int idir, int comp, amrex::Array4<amrex::Real const> const& q,
  amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& dxinv) noexcept
{
  using amrex::Real;
  const auto ivn = amrex::IntVect::TheDimensionVector(idir);

  return (fcw6_5 * q(iv + 2 * ivn, comp) + fcw6_4 * q(iv +     ivn, comp) + fcw6_3 * q(iv, comp) + 
          fcw6_2 * q(iv -     ivn, comp) + fcw6_1 * q(iv - 2 * ivn, comp)  + 
          fcw6_0 * q(iv - 3 * ivn, comp)) * dxinv[idir];
}

// Normal derivative (CELL CENTRED) with variable order using central differences ////
/// ********************************************************************** ////
template <int order>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real normal_diff_cc(
  amrex::IntVect iv, int idir, int comp, amrex::Array4<amrex::Real const> const& q,
  amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& dxinv) noexcept;
// cell centred
/// @brief 2nd order central difference face centred
template <>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real normal_diff_cc<2>(
  amrex::IntVect iv, int idir, int comp, amrex::Array4<amrex::Real const> const& q,
  amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& dxinv) noexcept
{
  const auto ivn = amrex::IntVect::TheDimensionVector(idir);
  return (ccw2_1*q(iv+ivn, comp) - ccw2_0*q(iv - ivn, comp)) * dxinv[idir];
}

/// @brief 4th order central difference face centred
template <>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real normal_diff_cc<4>(
  amrex::IntVect iv, int idir, int comp, amrex::Array4<amrex::Real const> const& q,
  amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& dxinv) noexcept
{
  using amrex::Real;
  const auto ivn = amrex::IntVect::TheDimensionVector(idir);

  return (ccw4_3 * q(iv + 2*ivn, comp) + ccw4_2 * q(iv + ivn, comp) +
          ccw4_1 * q(iv - ivn, comp) + ccw4_0 * q(iv - 2 * ivn, comp)) * dxinv[idir];
}

/// @brief 6th order central difference face centred
template <>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real normal_diff_cc<6>(
  amrex::IntVect iv, int idir, int comp, amrex::Array4<amrex::Real const> const& q,
  amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& dxinv) noexcept
{
  using amrex::Real;
  const auto ivn = amrex::IntVect::TheDimensionVector(idir);

  return (ccw6_5*q(iv + 3*ivn, comp) +  ccw6_4 * q(iv + 2*ivn, comp) + ccw6_3 * q(iv + ivn, comp) + 
          ccw6_2 * q(iv - ivn, comp) +  ccw6_1 * q(iv - 2*ivn, comp) + ccw6_0 * q(iv - 3 * ivn, comp)) 
          * dxinv[idir];
}

/// ********************************************************************** ////



//// Tangential derivative with variable order using central difference. ////

template <int order>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
tangent_diff(amrex::IntVect iv, int idir, int tdir, int comp,
             amrex::Array4<amrex::Real const> const& q,
             amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& dxinv) noexcept;

/// @brief 2nd order central difference. d/d(tdir) on face normal to idir.
template <>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
tangent_diff<2>(amrex::IntVect iv, int idir, int tdir, int comp,
                amrex::Array4<amrex::Real const> const& q,
                amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& dxinv) noexcept
{
  const auto ivn = amrex::IntVect::TheDimensionVector(idir);
  const auto ivt = amrex::IntVect::TheDimensionVector(tdir);

  // const auto ivm1 = iv - ivn;
  // const Real d2 = ccw2_1*q(iv + ivt, comp)     + ccw2_0*q(iv - ivt, comp);
  // const Real d1 = ccw2_1*q(ivm1   + ivt, comp) + ccw2_0*q(ivm1 - ivt, comp);
  // return( (d2*int_fcw2_1 + d1*int_fcw2_0) * dxinv[tdir]);

  return (q(iv + ivt, comp) - q(iv - ivt, comp) + q(iv - ivn + ivt, comp) -
          q(iv - ivn - ivt, comp)) * 0.25 * dxinv[tdir];
     

}

/// @brief 4th order central difference. d/d(tdir) on face normal to idir.
template <>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
tangent_diff<4>(amrex::IntVect iv, int idir, int tdir, int comp,
                amrex::Array4<amrex::Real const> const& q,
                amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& dxinv) noexcept
{
  using amrex::Real;
  const auto ivn = amrex::IntVect::TheDimensionVector(idir);  // for example i  - x
  const auto ivt = amrex::IntVect::TheDimensionVector(tdir);  //             j  - y

  // calculate cell-centred derivatives
  // dqdy (at i+1) using cell-centred derivatives(using j+2,j+1,j-1,j-2)
  const Real d3 = ccw4_3 * q(iv + ivn + 2 * ivt, comp) + ccw4_2 * q(iv + ivn + ivt, comp) +
            ccw4_1 * q(iv + ivn - ivt, comp)     + ccw4_0 * q(iv + ivn - 2 * ivt, comp);
  // dqdy (at i)          
  const Real d2 = ccw4_3 * q(iv + 2 * ivt, comp) + ccw4_2 * q(iv + ivt, comp) +
            ccw4_1  * q(iv - ivt, comp)    + ccw4_0 * q(iv - 2 * ivt, comp);
  // dqdy (at i -1)          
  const Real d1 = ccw4_3 * q(iv - ivn + 2 * ivt, comp) + ccw4_2 * q(iv - ivn + ivt, comp) +
            ccw4_1 * q(iv - ivn - ivt, comp)     + ccw4_0 * q(iv - ivn - 2 * ivt, comp);
  // dqdy (at i -2) 
  const Real d0 = ccw4_3 * q(iv - 2 * ivn + 2 * ivt, comp) + ccw4_2 * q(iv - 2 * ivn + ivt, comp) +
            ccw4_1 * q(iv - 2 * ivn - ivt, comp)     + ccw4_0 * q(iv - 2 * ivn - 2 * ivt, comp);

  // Interpolate derivatives using face-centred weights 
  return (d3 * int_fcw4_3 + d2 *int_fcw4_2 + d1*int_fcw4_1 + d0*int_fcw4_0) * dxinv[tdir];
}
//---------------------------------------------------------------------------------------
/// @brief 6th order central difference. d/d(tdir) on face normal to idir TODO
template <>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
tangent_diff<6>(amrex::IntVect iv, int idir, int tdir, int comp,
                amrex::Array4<amrex::Real const> const& q,
                amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& dxinv) noexcept
{
  using amrex::Real;
  const auto ivn = amrex::IntVect::TheDimensionVector(idir);  // for example i  - x
  const auto ivt = amrex::IntVect::TheDimensionVector(tdir);  //             j  - y

  const auto ip2 =  iv + 2*ivn; //i + 2
  const auto ip1 =  iv +   ivn; //i + 1
  const auto im1 =  iv -   ivn; //i - 1
  const auto im2 =  iv - 2*ivn; //i - 2
  const auto im3 =  iv - 3*ivn; //i - 3
  

  // calculate cell-centred derivatives

  // dqdy (at i+2) using cell-centred derivatives(using j+3,j+2,j+1,j-1,j-2,j-3)
  const Real d5 = ccw6_5 * q(ip2 + 3 * ivt, comp) + ccw6_4 * q(ip2 + 2 * ivt, comp) +
                  ccw6_3 * q(ip2 +     ivt, comp) + ccw6_2 * q(ip2 -     ivt, comp) +
                  ccw6_1 * q(ip2 - 2 * ivt, comp) + ccw6_0 * q(ip2 - 3 * ivt, comp);
  // dqdy (at i+1) using cell-centred derivatives
  const Real d4 = ccw6_5 * q(ip1 + 3 * ivt, comp) + ccw6_4 * q(ip1 + 2 * ivt, comp) +
                  ccw6_3 * q(ip1 +     ivt, comp) + ccw6_2 * q(ip1 -     ivt, comp) +
                  ccw6_1 * q(ip1 - 2 * ivt, comp) + ccw6_0 * q(ip1 - 3 * ivt, comp);
  // dqdy (at i)          
  const Real d3 = ccw6_5 * q(iv  + 3 * ivt, comp) + ccw6_4 * q(iv  + 2 * ivt, comp) +
                  ccw6_3 * q(iv  +     ivt, comp) + ccw6_2 * q(iv  -     ivt, comp) +
                  ccw6_1 * q(iv  - 2 * ivt, comp) + ccw6_0 * q(iv  - 3 * ivt, comp);
  // dqdy (at i -1)          
  const Real d2 = ccw6_5 * q(im1 + 3 * ivt, comp) + ccw6_4 * q(im1 + 2 * ivt, comp) +
                  ccw6_3 * q(im1 +     ivt, comp) + ccw6_2 * q(im1 -     ivt, comp) +
                  ccw6_1 * q(im1 - 2 * ivt, comp) + ccw6_0 * q(im1 - 3 * ivt, comp);
  // dqdy (at i -2) 
  const Real d1 = ccw6_5 * q(im2 + 3 * ivt, comp) + ccw6_4 * q(im2 + 2 * ivt, comp) +
                  ccw6_3 * q(im2 +     ivt, comp) + ccw6_2 * q(im2 -     ivt, comp) +
                  ccw6_1 * q(im2 - 2 * ivt, comp) + ccw6_0 * q(im2 - 3 * ivt, comp);

  // dqdy (at i -3)
  const Real d0 = ccw6_5 * q(im3 + 3 * ivt, comp) + ccw6_4 * q(im3 + 2 * ivt, comp) +
                  ccw6_3 * q(im3 +     ivt, comp) + ccw6_2 * q(im3 -     ivt, comp) +
                  ccw6_1 * q(im3 - 2 * ivt, comp) + ccw6_0 * q(im3 - 3 * ivt, comp);
          
  // Interpolate derivatives using face-centred weights 
  return (d5 * int_fcw6_5 + d4 * int_fcw6_4 + d3 * int_fcw6_3 + 
          d2 * int_fcw6_4 + d1 * int_fcw6_3 + d0 * int_fcw6_0) * dxinv[tdir];
}

////----------Interpolation with variable order using polynomial -------------- ////
///  interpolates in face between i and i-1   (i-1/2)
template <int order>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
interp(amrex::IntVect iv, int idir, int comp,
       amrex::Array4<amrex::Real const> const& q) noexcept;

template <>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
interp<2>(amrex::IntVect iv, int idir, int comp,
          amrex::Array4<amrex::Real const> const& q) noexcept
{
  const auto ivn = amrex::IntVect::TheDimensionVector(idir);
  return 0.5 * (q(iv, comp) + q(iv - ivn, comp));
}

template <>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
interp<4>(amrex::IntVect iv, int idir, int comp,
          amrex::Array4<amrex::Real const> const& q) noexcept
{
  const auto ivn = amrex::IntVect::TheDimensionVector(idir);
  return (int_fcw4_3 * q(iv + ivn, comp) + int_fcw4_2 * q(iv, comp) + int_fcw4_1 * q(iv - ivn, comp) +
          int_fcw4_0 * q(iv - 2 * ivn, comp));
}

template <>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
interp<6>(amrex::IntVect iv, int idir, int comp,
          amrex::Array4<amrex::Real const> const& q) noexcept
{
  const auto ivn = amrex::IntVect::TheDimensionVector(idir);
  return (int_fcw6_5 * q(iv + 2*ivn, comp) + int_fcw6_4 * q(iv + ivn,  comp)  + int_fcw6_3 * q(iv, comp) + 
          int_fcw6_2 * q(iv -ivn, comp)    + int_fcw6_1 * q(iv - 2*ivn, comp) + int_fcw6_0 * q(iv - 3 * ivn, comp));

}



#endif