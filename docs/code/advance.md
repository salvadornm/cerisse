# Advance the solution

This page explains how the ```advance.cpp``` work (located in ```src/tim```), 
that forms the core of the source code
and calls the appropiate subroutines. Following the general form:

$$
\mbox{RHS}(U) = \mbox{Euler} + \mbox{Diffusive} + \mbox{Source}
$$

The **RHS** is computed in the ```compute_rhs``` function, 
while the Runge-Kutta steps are in ```advance```function (that calls ```compute_rhs```).

## Advance Solution

This function advance the solution a full time step, given a time step.

```cpp
Real CNS::advance(Real time, Real dt, int /*iteration*/, int /*ncycle*/) {
```

## Compute RHS

```cpp
void CNS::compute_rhs(MultiFab& statemf, Real dt, 
          FluxRegister* fr_as_crse, FluxRegister* fr_as_fine) {
```

This function loop inside over all boxes (or *fabs*) in the current computing structure (either GPU or CPU)

```cpp
for (MFIter mfi(statemf, false); mfi.isValid(); ++mfi) {
  Array4<Real> const& state = statemf.array(mfi);
  .. 
```

Two short-lived ArrayBox are created, the primitive ```prims``` and a temporary flux (```fluxt```)

```cpp
  FArrayBox primf(bxg, cls_h.NPRIM, The_Async_Arena());
  FArrayBox tempf(bxg, cls_h.NCONS, The_Async_Arena());
  Array4<Real> const& fluxt = tempf.array();
  Array4<Real> const& prims= primf.array();
```

The primitive variables are computed from the state (which store the conservative variables)

```cpp
  cls_h.cons2prims(mfi, state, prims);
```


### Solid Markers

In case solid geometries are used, the solid markers  are computed either from **IBM** or **EBM** classes from their respective multifabs and current refinement levels.
The markers are then combined into one  short- lived boolean array ```geoMarkers```

### Convective and Diffusive fluxes

The following lines computes:

$$
\mbox{RHS}(U)  \leftarrow   \mbox{Euler} + \mbox{Diffusive} 
$$

```cpp
#ifdef AMREX_USE_GPIBM    
    prob_rhs.eflux_ibm(geom, mfi, prims, fluxt, state, cls_d, geoMarkers);
    prob_rhs.dflux_ibm(geom, mfi, prims, fluxt, state, cls_d, geoMarkers);
#else
    prob_rhs.eflux(geom, mfi, prims, fluxt, state, cls_d);
    prob_rhs.dflux(geom, mfi, prims, fluxt, state, cls_d);
#endif
```

These lines overwrite the ```state``` with the RHS value, and the fluxes themselves only depend on the primitive array.  The ```eflux``` and ```dflux```functions are defined  in the flux definition(see [fluxsolver](fluxsolver.md) ). The state must be set to 0 in the first call.

### Wall fluxes

In case of interal solid boundaries coded by **EB** the wall fluxes need to be computed. 
This is activated only if the box contains the boundary

```cpp
  FabType t = flag.getType(ebbox);
  if (FabType::singlevalued == t){
    EBM::eb.ebflux(geom,mfi, prims, fluxt,state, cls_d,level);
  }
```

The ```ebflux``` function is computed in the ```ebm.h`` class and updates state (the RHS) 
in the cells that contain a boundary. Particular cell treaments

### Source Terms

The following line calculates the source term, which may include chemical source terms or user-defined terms, such as external forces.

```cpp
  prob_rhs.src(mfi, prims, state, cls_d, dt);
```
$$
\mbox{RHS}(U)  \leftarrow  \mbox{Source}
$$


### Zero-th solid points

The last part of the subroutine, makes sure that if EB/IBM methods are used, the RHS is 
0 for solid points.

```cpp
  amrex::ParallelFor(bx, cls_h.NCONS,
    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
      state(i,j,k,n) = state(i,j,k,n)*(1 - int(geoMarkers(i,j,k,0)));
    });
```
