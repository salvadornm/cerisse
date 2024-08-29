# Understanding  a numerical solver

This page expalins the rationale between a new finite volume sovler. 
The Skew-symmetric will be used as an example

## Finite volume solver

A naive finite volume solver would be, in 1D 
$$
\left. \frac{\partial U}{\partial t} \right|_{i} = ( F_{i+1/2} -  F_{i-1/2} )
$$

Where  $$U$ = ( \rho, \rho u, \rho e_T) $ the vector of conservative variables and 
$$F = ( \rho u, \rho u^2 + P, \rho u e_T + u P  )$$ would be a flux function.

The main looop is

```cpp
    // loop over directions
    for (int dir = 0; dir < AMREX_SPACEDIM; dir++) 
      {
      GpuArray<int, 3> vdir = {int(dir == 0), int(dir == 1), int(dir == 2)};
       int Qdir =  cls_t::QRHO + dir + 1; 
      // compute interface fluxes at i-1/2, j-1/2, k-1/2
      ParallelFor(bxgnodal,
                  [=,*this] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    this->flux_dir(i, j, k,Qdir, vdir, cons, prims, lambda_max, flx, cls);
                  });
      // add flux derivative to rhs, i.e.  rhs[n] + = (fi[n] - fi+1[n])/dx
      ParallelFor(bx, cls_t::NCONS,
                  [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
                    rhs(i, j, k, n) +=
                        dxinv[dir] * (flx(i, j, k, n) - flx(i+vdir[0], j+vdir[1], k+vdir[2], n));
                  });
      }
```
where the interface fluxes are computed, note that flux(i+1) is the flux at interface (i+1/2) following standard finite volume notation. 
The two calls compute first the flux and the computes the difference to store in the rhs array that would be used to advance the solution within the Runge-Kutta step.

`vdir` is an integer array that depends on the direction  `vdir (idir=0)= [1,0,0]`
and  `vdir(idir=1)=[0,1,0]` and so on. 
`Qdir` is an integer that tells in which position in the `prims` array (primitive variables array) the velocity corrdinate that corresponds to this direction is stored. For example in direction 0 (that corresponds to x axis), the x-velocity is stored in position `Qdir=QRHO+1` (1 more position relative to the position of density), in direction 1 (y-axis) the y-velocity is stored in `Qdir=QRHO+2` (2 positions away from density).
Velocity is always follows the density in the vraibale array.

## Calculate Flux 

```cpp
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void flux_dir(
    int i, int j, int k, int Qdir,const GpuArray<int, 3>& vdir, const Array4<Real>& cons, const Array4<Real>& prims, const Array4<Real>& lambda_max, const Array4<Real>& flx,
    const cls_t* cls) const {
    
        // prepare arrays
    int il= i-order_skew*vdir[0]; int jl= j-order_skew*vdir[1]; int kl= k-order_skew*vdir[2];   
    for (int l = 0; l < order_skew; l++) {  
      il +=  vdir[0];jl +=  vdir[1];kl +=  vdir[2];
      rho =prims(il,jl,kl,cls_t::QRHO);
      V[l] = prims(il,jl,kl,Qdir); 
      U[l][cls_t::QRHO] = rho;
      U[l][cls_t::QU] = rho*prims(il,jl,kl,cls_t::QU);
      U[l][cls_t::QV] = rho*prims(il,jl,kl,cls_t::QV);
      U[l][cls_t::QW] = rho*prims(il,jl,kl,cls_t::QW);
      eint=cls->cv * prims(il, jl, kl, cls_t::QT);
      U[l][cls_t::NCONS-1] = rho*(eint + kin) + prims(il,jl,kl,cls_t::QPRES);
    }

    ....

    }
```    
In this example, the flux is a  skew-symmetric formulation, in a second order is

$$
F_{i-1/2} = \frac{1}{2} \left( U(i-1) + U(i) \right) \frac{1}{2} \left( V(i-1) + V(i) \right)
$$


`il`,`jl` and `kl` are `i-1`,`j-1` and `k-1` but using the direction, so if the function is called within `idir=0` (x-direction), `il=i-1`,`jl=j` and `kl=k`, similarly with 
`idir=1` (y-direction), `il=i1`,`jl=j-1` and `kl=k`.
The `prepare arrays` section will build the conservative variable arrays 
$$ U_{i-order/2} , ... , U_{i+order/2} $$ 




